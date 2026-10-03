#include "photon_flux.hpp"
#include "def.hpp"
#include <memory>
#include <cmath>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_bessel.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

static double Kn_flux(int nu, double x)
{
    gsl_sf_result result;
    int status = gsl_sf_bessel_Kn_e(nu, x, &result);
    if (status != GSL_SUCCESS) return 0.0;
    return result.val;
}

static double get_GammaAA(double b) { return gamma_aa()(b); }

double flux_density(double qp, double r, void* p)
{
    parameters* par = (parameters*)p;
    double omega = qp / sqrt(2.0);
    double z_gamma = omega * 2.0 / par->ss;
    double eta     = z_gamma * par->mn * r;
    if (eta > 50.0) return 0.0;
    double pref = (par->alpha * par->Z * par->Z) / (M_PI * M_PI);
    double K1   = Kn_flux(1, eta);
    return (pref / omega) * (eta*eta / (r*r)) * (K1*K1);
}

namespace {
    std::unique_ptr<TASpline>     ta_spline;
    std::unique_ptr<GammaAA>      gamma_instance;
    std::unique_ptr<WSFormFactor> ws_form_factor;
    std::unique_ptr<WSFluxTable>  ws_flux_table;
    double ws_RA = 0.0;       // GeV^-1, set by init_ws_form_factor()
    double ws_r_switch = 0.0; // GeV^-1: below it use the WS table, above it the point-like flux (default 2R_PL)
    const double WS_QMAX = 5.0;   // GeV, see init_ws_form_factor()
}

// Woods-Saxon photon flux. For r >= 2R_PL = 14.2 fm it is the same as the
// point-like flux_density(), as in the appendix of arXiv:2404.09731.
double flux_density_WS(double z_gamma, double r, void* p)
{
    parameters* par = (parameters*)p;
    if (!ws_flux_table)
        throw std::runtime_error("flux_density_WS() called before init_ws_form_factor()");
    if (r >= ws_r_switch) {
        double omega = z_gamma * par->ss / 2.0;
        return flux_density(omega * std::sqrt(2.0), r, p);
    }

    // The factor 2/ss makes this equal to flux_density() when F_A = 1.
    double pref = (2.0 * par->alpha * par->Z * par->Z) / (M_PI * M_PI * z_gamma * par->ss);
    return pref * ws_flux_table->I2(z_gamma, r);
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
    double z_gamma = omega * 2.0 / par->ss;
    // the flux is taken at the centre of the target nucleus, r = b
    double r     = b;
    double flux  = (par->flux_model == "WS") ? flux_density_WS(z_gamma, r, p)
                                              : flux_density(qp, r, p);

    if (par->target == "pA") {
        // pA: no EMD factor. "PL" is zero below bmin. The others use Gamma_pA(b).
        double gamma = (par->flux_model == "PL") ? ((b >= par->bmin) ? 1.0 : 0.0)
                                                  : GammaPA(b);
        if (gamma <= 0.0) return 0.0;
        return 2.0 * M_PI * b * flux * gamma;
    }

    // only the point-like channel is zero below bmin (Gamma = 1 above it)
    double gamma = (channel == "PL(AnAn)") ? ((b >= par->bmin) ? 1.0 : 0.0)
                                           : ((b < 150.0) ? get_GammaAA(b) : 1.0);
    if (gamma <= 0.0) return 0.0;   // this also avoids inf*0 for the point-like flux at b -> 0
    double P_emd      = par->S / (b * b);
    double P_no_emd   = std::exp(-P_emd);
    double emd_factor = 1.0;
    if      (channel == "An0n") emd_factor = P_no_emd;
    else if (channel == "Xn0n") emd_factor = P_no_emd * (1.0 - P_no_emd);
    else if (channel == "0n0n") emd_factor = P_no_emd * P_no_emd;   // no nucleus breaks up
    return 2.0 * M_PI * b * flux * gamma * emd_factor;
}

namespace {
    std::string gamma_aa_source;   // file (and sigma_NN) of the loaded Gamma_AA table
}

const std::string& gamma_aa_info() { return gamma_aa_source; }

void load_data_and_initialize(const std::string& default_filename)
{
    // GAMMA_AA_FILE: use another Gamma_AA table.
    const char* env = std::getenv("GAMMA_AA_FILE");
    const std::string filename = env ? env : default_filename;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + filename);
    }

    std::vector<double> b_values;
    std::vector<double> T_values;
    std::string line;
    gamma_aa_source = filename;

    // columns: b [GeV^-1]  Gamma_AA(b)
    while (std::getline(file, line)) {
        // the header line "# sigma_NN : 92 mb" says which sigma_NN the table was made with
        if (line.rfind("# sigma_NN", 0) == 0) {
            gamma_aa_source += " (sigma_NN = " + line.substr(line.find(':') + 1) + ")";
            gamma_aa_source.erase(gamma_aa_source.find("=  ") + 2, 1);
        }
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        double b, T;
        if (iss >> b >> T) {
            b_values.push_back(b);
            T_values.push_back(T);
        }
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
    // radius where we change from WS to PL, 2R_PL = 14.2 fm (change it with WS_RSWITCH_FM)
    ws_r_switch = (std::getenv("WS_RSWITCH_FM") ? std::atof(std::getenv("WS_RSWITCH_FM")) : 14.2) / hbarc;
    double a_GeV = a_fm / hbarc;
    ws_form_factor.reset(new WSFormFactor(ws_RA, a_GeV, WS_QMAX));
    // The r grid stops at that radius (above it we use the point-like flux).
    // Grid sizes can be changed for tests with WS_NZ, WS_NB.
    int nz_gamma = std::getenv("WS_NZ") ? std::atoi(std::getenv("WS_NZ")) : 160;
    int nr = std::getenv("WS_NB") ? std::atoi(std::getenv("WS_NB")) : 260;
    ws_flux_table.reset(new WSFluxTable(*ws_form_factor, mn, WS_QMAX,
                                         1e-6, 1.0, nz_gamma,      // z_gamma grid
                                         1.0, ws_r_switch, nr));   // r grid
}

namespace {
    std::unique_ptr<WSThickness>   eff_TB;
    std::unique_ptr<EffFluxRadial> eff_H;
    gsl_spline* eff_flux_spline = nullptr;
    gsl_interp_accel* eff_flux_acc = nullptr;
    double eff_lz_gamma_min = 0.0, eff_lz_gamma_max = 0.0;
    const double EFF_R_SWITCH = 800.0;   // GeV^-1. With 300 the result was off by about 1%
}

void init_effective_flux(const std::string& channel, void* p, double RA_fm, double a_fm, double B_mass)
{
    parameters* par = (parameters*)p;
    init_ws_form_factor(RA_fm, a_fm, par->mn);   // the r integral below needs the flux down to r=0

    const double hbarc = 0.197327;
    double RA = RA_fm / hbarc, a = a_fm / hbarc;
    double s_max = RA + 12.0 * a;

    eff_TB.reset(new WSThickness(RA, a, B_mass, s_max));
    eff_H.reset(new EffFluxRadial(*eff_TB, channel, par->S, B_mass, EFF_R_SWITCH));

    // f_eff(z_gamma) = (2*pi/B) \int_0^{r_max(z_gamma)} r dr f_WS(z_gamma,r) H_channel(r), r_max = 60/(z_gamma*mn)
    int nz_gamma = 100;
    double z_gamma_min = 1e-6, z_gamma_max = 1.0;
    eff_lz_gamma_min = std::log(z_gamma_min); eff_lz_gamma_max = std::log(z_gamma_max);
    std::vector<double> lz_gamma(nz_gamma), logf(nz_gamma);
    for (int i = 0; i < nz_gamma; i++) {
        double lz_gamma_i = eff_lz_gamma_min + (eff_lz_gamma_max - eff_lz_gamma_min) * i / (nz_gamma - 1.0);
        lz_gamma[i] = lz_gamma_i;
        double z_gamma = std::exp(lz_gamma_i);
        double r_max = 60.0 / (z_gamma * par->mn);
        double r0 = 1e-3;
        int Nr = 2000;
        double aa = std::log(r0), cc = std::log(r_max), hh = (cc - aa) / Nr, sum = 0.0;
        for (int k = 0; k <= Nr; k++) {
            double r = std::exp(aa + k * hh);
            double flux = flux_density_WS(z_gamma, r, p);
            double Hval = (*eff_H)(r);
            double w = (k == 0 || k == Nr) ? 1.0 : (k % 2 ? 4.0 : 2.0);
            sum += w * (r * flux * Hval) * r;   // extra r: dr = r dlnr
        }
        double r_integral = sum * hh / 3.0;
        double f_eff = (2.0 * M_PI / B_mass) * r_integral;
        // NaN or 0 (at z_gamma = 1, where the flux vanishes): continue the line of the two points
        // before, since a jump to a tiny value makes the spline oscillate up to z_gamma ~ 0.5
        if (std::isfinite(f_eff) && f_eff > 0.0) logf[i] = std::log(f_eff);
        else logf[i] = (i >= 2) ? 2.0 * logf[i-1] - logf[i-2] : std::log(1e-300);
    }
    if (eff_flux_spline) { gsl_spline_free(eff_flux_spline); gsl_interp_accel_free(eff_flux_acc); }
    eff_flux_spline = gsl_spline_alloc(gsl_interp_cspline, nz_gamma);
    eff_flux_acc = gsl_interp_accel_alloc();
    gsl_spline_init(eff_flux_spline, lz_gamma.data(), logf.data(), nz_gamma);
}

namespace {
    std::vector<double> sl_x, sl_logf;   // y^(1/4) points and log dN/dy (Starlight table)
    double sl_A = 0.0, sl_logc = 0.0;    // used below the first table point, see init_starlight_flux()
}

void init_starlight_flux(const std::string& channel, const std::string& filename)
{
    int col;
    if      (channel == "AnAn") col = 3;
    else if (channel == "An0n") col = 4;
    else throw std::runtime_error("FLUX_MODEL=STARLIGHT: the table only has AnAn and An0n, not " + channel);

    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Could not open file " + filename);
    sl_x.clear(); sl_logf.clear();
    std::string line;
    // columns: j  y  y^(1/4)  log_AnAn  log_An0n
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        double c[5];
        if (!(iss >> c[0] >> c[1] >> c[2] >> c[3] >> c[4])) continue;
        sl_x.push_back(c[2]);
        sl_logf.push_back(c[col]);
    }
    if (sl_x.size() < 2) throw std::runtime_error("Starlight table " + filename + " has no data");

    // Below the first table point (y = 1e-4) use y*f(y) = A ln(c/y),
    // fixed by the first two table points.
    double y0 = std::pow(sl_x[0], 4), y1 = std::pow(sl_x[1], 4);
    double g0 = y0 * std::exp(sl_logf[0]), g1 = y1 * std::exp(sl_logf[1]);
    sl_A    = (g0 - g1) / std::log(y1 / y0);
    sl_logc = g0 / sl_A + std::log(y0);
}

// dN/dy from the table, interpolated in y^(1/4) with the formula of
// Eqs. 34-37 of arXiv:2404.09731.
static double starlight_dNdy(double y)
{
    if (y >= 1.0) return 0.0;
    double y0 = std::pow(sl_x.front(), 4);
    if (y < y0) return sl_A * (sl_logc - std::log(y)) / y;
    double x = std::pow(y, 0.25), num = 0.0, den = 0.0;
    int n = (int)sl_x.size() - 1;
    for (int j = 0; j <= n; j++) {
        double d = x - sl_x[j];
        if (std::fabs(d) < 1e-14) return std::exp(sl_logf[j]);
        double w = ((j % 2) ? -1.0 : 1.0) * ((j == 0 || j == n) ? 0.5 : 1.0) / d;
        num += w * sl_logf[j];
        den += w;
    }
    return std::exp(num / den);
}

double effective_photon_flux(double qp, void* p)
{
    parameters* par = (parameters*)p;
    if (par->flux_model == "STARLIGHT") {
        if (sl_x.empty())
            throw std::runtime_error("effective_photon_flux() called before init_starlight_flux()");
        double z_gamma = (qp / std::sqrt(2.0)) * 2.0 / par->ss;
        return starlight_dNdy(z_gamma) * 2.0 / par->ss;   // dN/domega, same units as f_eff
    }
    if (!eff_flux_spline)
        throw std::runtime_error("effective_photon_flux() called before init_effective_flux()");
    double omega = qp / std::sqrt(2.0);
    double z_gamma = omega * 2.0 / par->ss;
    double lz_gamma = std::log(z_gamma);
    if (lz_gamma < eff_lz_gamma_min) lz_gamma = eff_lz_gamma_min; else if (lz_gamma > eff_lz_gamma_max) lz_gamma = eff_lz_gamma_max;
    return std::exp(gsl_spline_eval(eff_flux_spline, lz_gamma, eff_flux_acc));
}