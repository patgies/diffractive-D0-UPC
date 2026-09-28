#include "photon_flux.hpp"
#include "def.hpp"
#include <memory>
#include <cmath>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

static double Kn_flux(int nu, double x)
{
    gsl_sf_result result;
    int status = gsl_sf_bessel_Kn_e(nu, x, &result);
    if (status != GSL_SUCCESS) return 0.0;
    return result.val;
}

static double get_GammaAA(double b) { return gamma_aa()(b); }

double flux_density(double qp, double b, void* p)
{
    parameters* par = (parameters*)p;
    double omega = qp / sqrt(2.0);
    double z     = omega * 2.0 / par->ss;
    double eta   = z * par->mn * b;
    if (eta > 50.0) return 0.0;
    double pref = (par->alpha * par->Z * par->Z) / (M_PI * M_PI);
    double K1   = Kn_flux(1, eta);
    return (pref / omega) * (eta*eta / (b*b)) * (K1*K1);
}

namespace {
    std::unique_ptr<TASpline>     ta_spline;
    std::unique_ptr<GammaAA>      gamma_instance;
    std::unique_ptr<WSFormFactor> ws_form_factor;
    std::unique_ptr<WSFluxTable>  ws_flux_table;
    double ws_RA = 0.0;       // GeV^-1, set by init_ws_form_factor()
    const double WS_QMAX = 5.0;   // GeV -- see init_ws_form_factor()
}

// Woods-Saxon bare photon flux (Vidovic-Greiner-Best-Soff / Krauss-Greiner-Soff;
// see photon_flux.hpp for WSFormFactor/WSFluxTable and the derivation
// references). Falls back to the exact closed-form flux_density() for
// b >= R_A: past the nuclear radius, the WS bare flux is IDENTICAL to the
// point-like one (verified numerically to double precision; see
// init_ws_form_factor()), so this is only ever different from flux_density()
// for b < R_A -- and since par->bmin is chosen as ~2*R_A everywhere in this
// code (to exclude nuclear overlap), FLUX_MODEL=WS and FLUX_MODEL=PL give
// IDENTICAL cross sections in the current pipeline (b < R_A is physically
// excluded there anyway). The source nucleus's own charge form factor is
// therefore NOT what drives the PL-vs-"Starlight" discrepancy seen in
// flux_comparison.pdf -- that must come from the target-nucleus geometric
// convolution (the "effective flux", Eq. 4 of arXiv:2404.09731), which this
// function does not implement.
double flux_density_WS(double z, double b, void* p)
{
    parameters* par = (parameters*)p;
    if (!ws_flux_table)
        throw std::runtime_error("flux_density_WS() called before init_ws_form_factor()");
    if (b >= ws_RA) {
        double omega = z * par->ss / 2.0;
        return flux_density(omega * std::sqrt(2.0), b, p);
    }

    // The literature Eq. 16-style prefactor is alpha*Z^2/(pi^2 z), but this
    // code's own flux_density() (validated against the tabulated PL flux)
    // uses a normalization convention that differs from the bare literature
    // formula by a factor of ss/2 (equivalently, per unit omega rather than
    // per unit z -- confirmed by matching this formula's own F_A=1 limit,
    // z*mN*K1(eta), against flux_density()'s eta^2/b^2/omega form). Include
    // the matching 2/ss here so flux_density_WS() reduces to flux_density()
    // exactly when F_A=1, i.e., so "WS" and "PL" only differ by the nuclear
    // form factor, nothing else.
    double pref = (2.0 * par->alpha * par->Z * par->Z) / (M_PI * M_PI * z * par->ss);
    return pref * ws_flux_table->I2(z, b);
}

namespace {
    std::unique_ptr<WSThickness> pA_TA;
    double pA_sigma_NN = 0.0;   // GeV^-2
}

void init_pA_flux(double sigma_NN_mb, double RA_fm, double a_fm, double B_mass)
{
    const double hbarc = 0.197327;
    double RA = RA_fm / hbarc, a = a_fm / hbarc;
    double s_max = RA + 12.0 * a;
    pA_TA.reset(new WSThickness(RA, a, B_mass, s_max));
    // 1 mb = 0.1 fm^2 = 0.1/hbarc^2 GeV^-2.
    pA_sigma_NN = sigma_NN_mb * 0.1 / (hbarc * hbarc);
}

double GammaPA(double b)
{
    if (!pA_TA)
        throw std::runtime_error("GammaPA() called before init_pA_flux()");
    return std::exp(-pA_sigma_NN * (*pA_TA)(b));
}

double photon_flux(double b, double qp, void* p)
{
    parameters* par = (parameters*)p;
    const string& channel = par->channel;
    double omega = qp / sqrt(2.0);
    double z     = omega * 2.0 / par->ss;
    double flux  = (par->flux_model == "WS") ? flux_density_WS(z, b, p)
                                              : flux_density(qp, b, p);

    if (par->target == "pA") {
        // Sec. 6.1 of arXiv:2606.05469: no EMD factor at all for pA. "WS"
        // (default, realistic) uses the smooth Gamma_pA(b)=exp(-sigma_NN*T_A(b))
        // survival factor; "PL" is the sharp point-like comparison cutoff at
        // par->bmin (set to 1.1*R_A by whoever configures target=="pA").
        double gamma = (par->flux_model == "PL") ? ((b >= par->bmin) ? 1.0 : 0.0)
                                                  : GammaPA(b);
        return 2.0 * M_PI * b * flux * gamma;
    }

    double gamma = (channel == "PL(AnAn)") ? ((b >= 14.0/0.197327) ? 1.0 : 0.0)
                                           : ((b < 150.0) ? get_GammaAA(b) : 1.0);
    double P_emd      = par->S / (b * b);
    double P_no_emd   = std::exp(-P_emd);
    double emd_factor = 1.0;
    if      (channel == "An0n") emd_factor = P_no_emd;
    else if (channel == "Xn0n") emd_factor = P_no_emd * (1.0 - P_no_emd);
    return 2.0 * M_PI * b * flux * gamma * emd_factor;
}

void load_data_and_initialize(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file!");
    }

    std::vector<double> b_values;
    std::vector<double> T_values;
    double b, T;

    while (file >> b >> T) {
        b_values.push_back(b);
        T_values.push_back(T);
    }


    ta_spline.reset(new TASpline(b_values, T_values));
    gamma_instance.reset(new GammaAA(*ta_spline));
}

GammaAA& gamma_aa()
{
    if (!gamma_instance) {
        throw std::runtime_error(
            "gamma_aa() called before load_data_and_initialize()");
    }
    return *gamma_instance;
}

void init_ws_form_factor(double RA_fm, double a_fm, double mn)
{
    const double hbarc = 0.197327;
    ws_RA = RA_fm / hbarc;
    double a_GeV = a_fm / hbarc;
    ws_form_factor.reset(new WSFormFactor(ws_RA, a_GeV, WS_QMAX));
    // Grid: z covers everything the flux is ever evaluated at (see
    // scan_flux.cpp/main.cpp). b only needs to reach R_A: for b >= R_A (i.e.
    // outside the emitting nucleus), the WS bare flux is IDENTICAL to the
    // point-like one -- confirmed numerically to ~1e-4 already at b=R_A and
    // to double precision by b=1.5*R_A -- the EPA analogue of Newton's shell
    // theorem (the field outside a spherically symmetric charge distribution
    // doesn't depend on its internal structure). flux_density_WS() falls
    // back to the exact closed-form flux_density() there, both because it's
    // exact and because the k_perp integral becomes numerically unstable
    // (catastrophic cancellation between the J1 oscillation and the nuclear
    // form factor's diffraction zero) once it's just reproducing that same
    // number from scratch.
    ws_flux_table.reset(new WSFluxTable(*ws_form_factor, mn, WS_QMAX,
                                         1e-6, 1.0, 80,      // z grid
                                         1.0, ws_RA, 60));   // b grid
}

namespace {
    std::unique_ptr<WSThickness>   eff_TB;
    std::unique_ptr<EffFluxRadial> eff_H;
    gsl_spline* eff_flux_spline = nullptr;
    gsl_interp_accel* eff_flux_acc = nullptr;
    double eff_lz_min = 0.0, eff_lz_max = 0.0;
    const double EFF_R_SWITCH = 800.0;   // GeV^-1 -- see EffFluxRadial's comment: 300 left a
                                          // percent-level mismatch for the EMD-bearing channels
}

void init_effective_flux(const std::string& channel, void* p, double RA_fm, double a_fm, double B_mass)
{
    parameters* par = (parameters*)p;
    init_ws_form_factor(RA_fm, a_fm, par->mn);   // the r-integral below needs the bare flux down to r=0

    const double hbarc = 0.197327;
    double RA = RA_fm / hbarc, a = a_fm / hbarc;
    double s_max = RA + 12.0 * a;

    eff_TB.reset(new WSThickness(RA, a, B_mass, s_max));
    eff_H.reset(new EffFluxRadial(*eff_TB, channel, par->S, par->bmin, B_mass, EFF_R_SWITCH));

    // f_eff(z) = (2*pi/B) \int_0^{r_max(z)} r dr f_gamma/A^WS(z,r) H_channel(r),
    // same r_max(z)=60/(z*mn) adaptive cutoff as scan_flux.cpp/main.cpp.
    int nz = 100;
    double z_min = 1e-6, z_max = 1.0;
    eff_lz_min = std::log(z_min); eff_lz_max = std::log(z_max);
    std::vector<double> lz(nz), logf(nz);
    for (int i = 0; i < nz; i++) {
        double lzz = eff_lz_min + (eff_lz_max - eff_lz_min) * i / (nz - 1.0);
        lz[i] = lzz;
        double z = std::exp(lzz);
        double r_max = 60.0 / (z * par->mn);
        double r0 = 1e-3;
        int Nr = 2000;
        double aa = std::log(r0), cc = std::log(r_max), hh = (cc - aa) / Nr, sum = 0.0;
        for (int k = 0; k <= Nr; k++) {
            double r = std::exp(aa + k * hh);
            double flux = flux_density_WS(z, r, p);
            double Hval = (*eff_H)(r);
            double w = (k == 0 || k == Nr) ? 1.0 : (k % 2 ? 4.0 : 2.0);
            sum += w * (r * flux * Hval) * r;   // extra r: dr = r dlnr
        }
        double r_integral = sum * hh / 3.0;
        double f_eff = (2.0 * M_PI / B_mass) * r_integral;
        logf[i] = std::log(std::max(f_eff, 1e-300));
    }
    if (eff_flux_spline) { gsl_spline_free(eff_flux_spline); gsl_interp_accel_free(eff_flux_acc); }
    eff_flux_spline = gsl_spline_alloc(gsl_interp_cspline, nz);
    eff_flux_acc = gsl_interp_accel_alloc();
    gsl_spline_init(eff_flux_spline, lz.data(), logf.data(), nz);
}

double effective_photon_flux(double qp, void* p)
{
    parameters* par = (parameters*)p;
    if (!eff_flux_spline)
        throw std::runtime_error("effective_photon_flux() called before init_effective_flux()");
    double omega = qp / std::sqrt(2.0);
    double z     = omega * 2.0 / par->ss;
    double lz = std::log(z);
    if (lz < eff_lz_min) lz = eff_lz_min; else if (lz > eff_lz_max) lz = eff_lz_max;
    return std::exp(gsl_spline_eval(eff_flux_spline, lz, eff_flux_acc));
}