#ifndef PHOTON_FLUX_HPP
#define PHOTON_FLUX_HPP
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>
#include <stdexcept>
#include <fstream>


// Woods-Saxon (Fermi-model) nuclear charge form factor F_A(q), the exact
// closed form of Maximon & Schrack, "The form factor of the Fermi model
// spatial distribution," J. Res. Natl. Bur. Stand. B 70, 85 (1966) --
// referenced in the technical appendix of arXiv:2404.09731 as the method
// used for their own tabulated flux (in place of a brute-force numerical
// Fourier transform, which is what this class used to do; the switch away
// from that is what this comment is about). For
//   rho_A(r) = rho0/(1+exp((r-c)/a))   (c = R_A, the half-density radius),
// their Eq. (20) gives (their Eq. (9) substitution qa=beta is not needed
// here -- this is the final, re-substituted result):
//   F(q) = [4 pi^2 rho0 a^3 / ((qa)^2 sinh^2(pi q a))]
//            * [pi q a cosh(pi q a) sin(qc) - qc cos(qc) sinh(pi q a)]
//          + 8 pi rho0 a^3 sum_{n=1}^inf (-1)^{n-1} n e^{-n c/a} / (n^2+(qa)^2)^2 ,
// with rho0 fixed by their Eq. (21) so that F(0)=1 exactly (no separate
// normalization step needed, unlike the old numerical version). The series
// converges extremely fast for realistic c/a (~12 for lead): a handful of
// terms is already far below double precision, so unlike the old version
// this needs no grid or spline at all -- it's evaluated in closed form
// every time, which is both exact (up to truncating the series) and cheap.
class WSFormFactor {
private:
    double c_, a_, q_max_, rho0_;
    static const int N_TERMS = 10;

    static double series_cube(double c, double a)   // sum (-1)^(n-1) e^{-nc/a}/n^3, for rho0
    {
        double sum = 0.0, sign = 1.0;
        for (int n = 1; n <= N_TERMS; n++) {
            sum += sign * std::exp(-n * c / a) / (n * n * n);
            sign = -sign;
        }
        return sum;
    }

    static double series_F(double q, double c, double a)   // sum (-1)^(n-1) n e^{-nc/a}/(n^2+(qa)^2)^2
    {
        double qa2 = (q * a) * (q * a);
        double sum = 0.0, sign = 1.0;
        for (int n = 1; n <= N_TERMS; n++) {
            double denom = (n * n + qa2);
            sum += sign * n * std::exp(-n * c / a) / (denom * denom);
            sign = -sign;
        }
        return sum;
    }

public:
    // RA, a in the same length units as the b/r used elsewhere in this code
    // (GeV^-1 here -- convert from fm with hbarc=0.197327 before calling).
    WSFormFactor(double RA, double a, double q_max) : c_(RA), a_(a), q_max_(q_max)
    {
        // Eq. (21): rho0 = { (4 pi c/3)[(pi a)^2+c^2] + 8 pi a^3 * series_cube(c,a) }^-1
        double bracket = (4.0 * M_PI * c_ / 3.0) * (M_PI * a_ * M_PI * a_ + c_ * c_)
                        + 8.0 * M_PI * a_ * a_ * a_ * series_cube(c_, a_);
        rho0_ = 1.0 / bracket;
    }

    WSFormFactor(const WSFormFactor&) = delete;
    WSFormFactor& operator=(const WSFormFactor&) = delete;

    double operator()(double q) const
    {
        if (q >= q_max_) return 0.0;      // fully decayed well before q_max by construction
        if (q < 1e-6) return 1.0;         // F(0)=1; removable 0/0 in the elementary term below
        double qa  = q * a_;
        double pqa = M_PI * qa;
        double sh = std::sinh(pqa), ch = std::cosh(pqa);
        double elementary = (4.0 * M_PI * M_PI * rho0_ * a_ * a_ * a_) / (qa * qa * sh * sh)
                           * (pqa * ch * std::sin(q * c_) - q * c_ * std::cos(q * c_) * sh);
        double series = 8.0 * M_PI * rho0_ * a_ * a_ * a_ * series_F(q, c_, a_);
        return elementary + series;
    }
};

#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_sf_result.h>

// Precomputed table of the Woods-Saxon bare-flux k_perp integral,
//   I(z,b) = \int_0^{q_max} dk_perp k_perp^2 F_A(sqrt(k_perp^2+z^2 mN^2))
//                            / (k_perp^2+z^2 mN^2) * J_1(k_perp b) ,
// so that f_{gamma/A}^WS(z,b) = (alpha_em Z^2 / (pi^2 z)) * I(z,b)^2 (Eq. 16
// of arXiv:2606.05469 / Eq. 11 of arXiv:2404.09731). J_1(k_perp b) oscillates
// too fast to integrate live at every Monte Carlo point once b is more than
// a few tens of GeV^-1, so I(z,b) is tabulated once (log z, log b grid) and
// looked up via a 2D spline -- same reason Petja Paakkinen's own tables
// (inputs/photon_flux/) are precomputed rather than evaluated on the fly.
// Only ever needed for b <= ~5*R_A (flux_density_WS() falls back to the
// point-like flux beyond that), which keeps the b-range -- and thus the
// number of J_1 oscillations per grid point -- manageable.
class WSFluxTable {
private:
    gsl_spline2d* spline;
    gsl_interp_accel* z_acc;
    gsl_interp_accel* b_acc;
    double lz_min, lz_max, lb_min, lb_max;

    struct IntegrandParams { double z, b, mn2; const WSFormFactor* FA; };

    static double integrand(double kt, void* p)
    {
        IntegrandParams* par = (IntegrandParams*)p;
        double t = kt * kt + par->z * par->z * par->mn2;
        gsl_sf_result res;
        int status = gsl_sf_bessel_Jn_e(1, kt * par->b, &res);
        double J1 = (status == GSL_SUCCESS) ? res.val : 0.0;
        return kt * kt * (*par->FA)(std::sqrt(t)) / t * J1;
    }

    static double raw_I(double z, double b, double mn, const WSFormFactor& FA, double q_max)
    {
        IntegrandParams par{z, b, mn * mn, &FA};
        gsl_function F; F.function = &integrand; F.params = &par;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(2000);
        double result, err;
        gsl_integration_qag(&F, 0.0, q_max, 0, 1e-5, 2000, GSL_INTEG_GAUSS61, w, &result, &err);
        gsl_integration_workspace_free(w);
        return result;
    }

public:
    WSFluxTable(const WSFormFactor& FA, double mn, double q_max,
                double z_min, double z_max, int nz,
                double b_min, double b_max, int nb)
    {
        lz_min = std::log(z_min); lz_max = std::log(z_max);
        lb_min = std::log(b_min); lb_max = std::log(b_max);

        std::vector<double> lz(nz), lb(nb);
        for (int i = 0; i < nz; i++) lz[i] = lz_min + (lz_max - lz_min) * i / (nz - 1.0);
        for (int j = 0; j < nb; j++) lb[j] = lb_min + (lb_max - lb_min) * j / (nb - 1.0);

        spline = gsl_spline2d_alloc(gsl_interp2d_bicubic, nz, nb);
        std::vector<double> grid_z(nz * nb);   // gsl_spline2d's own flattened layout
        for (int i = 0; i < nz; i++) {
            double z = std::exp(lz[i]);
            for (int j = 0; j < nb; j++) {
                double b = std::exp(lb[j]);
                double I = raw_I(z, b, mn, FA, q_max);
                // isfinite: a failed k_perp integral (seen near z=1) must not poison the spline
                double logI2 = std::isfinite(I) ? std::log(std::max(I * I, 1e-300)) : std::log(1e-300);
                gsl_spline2d_set(spline, grid_z.data(), i, j, logI2);
            }
        }
        gsl_spline2d_init(spline, lz.data(), lb.data(), grid_z.data(), nz, nb);
        z_acc = gsl_interp_accel_alloc();
        b_acc = gsl_interp_accel_alloc();
    }

    WSFluxTable(const WSFluxTable&) = delete;
    WSFluxTable& operator=(const WSFluxTable&) = delete;

    // Returns I(z,b)^2 (not yet multiplied by the alpha_em Z^2/(pi^2 z) prefactor).
    double I2(double z, double b) const
    {
        double lz_ = std::log(z), lb_ = std::log(b);
        if (lz_ < lz_min) lz_ = lz_min; else if (lz_ > lz_max) lz_ = lz_max;
        if (lb_ < lb_min) lb_ = lb_min; else if (lb_ > lb_max) lb_ = lb_max;
        return std::exp(gsl_spline2d_eval(spline, lz_, lb_, z_acc, b_acc));
    }

    ~WSFluxTable()
    {
        gsl_spline2d_free(spline);
        gsl_interp_accel_free(z_acc);
        gsl_interp_accel_free(b_acc);
    }
};


class TASpline {
private:
    gsl_spline* spline;
    gsl_interp_accel* acc;
    double s_min, s_max;
 
public:
    TASpline(const std::vector<double>& s,
             const std::vector<double>& T)
    {
        spline = gsl_spline_alloc(gsl_interp_cspline, s.size());
        acc = gsl_interp_accel_alloc();
        gsl_spline_init(spline, s.data(), T.data(), s.size());
 
        s_min = s.front();
        s_max = s.back();
    }
 
    
    TASpline(const TASpline&) = delete;
    TASpline& operator=(const TASpline&) = delete;
 
    double operator()(double s)
    {
        if (s < s_min || s > s_max)
            return 1.0;
 
        return gsl_spline_eval(spline, s, acc);
    }
 
    ~TASpline()
    {
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
    }
};
 

class GammaAA {
private:
    TASpline& TA;
    double b_min;
    double b_max;
 
public:
    GammaAA(TASpline& ta)
        : TA(ta), b_min(0.0), b_max(1e6)
    {}
 
    double operator()(double b)
    {
        if (b < b_min)
            throw std::out_of_range("b too small");
 
        if (b > b_max)
            return 1.0;
 
        // clamp: the cubic spline undershoots to ~-1e-32 where Gamma -> 0
        return std::min(1.0, std::max(0.0, TA(b)));
    }
};

GammaAA& gamma_aa();   // forward-declared here so EffFluxRadial below can use it

// Nuclear thickness function T_B(s): line-of-sight integral of the
// Woods-Saxon density, normalized so that 2*pi*\int_0^infty T_B(s) s ds =
// B_mass (number of nucleons). Needed for the "effective flux" geometric
// convolution (EffFluxRadial below) -- Eq. 4 of arXiv:2404.09731.
class WSThickness {
private:
    gsl_spline* spline;
    gsl_interp_accel* acc;
    double s_max_;

public:
    WSThickness(double RA, double a, double B_mass, double s_max, int n = 150)
        : s_max_(s_max)
    {
        double zmax = RA + 12.0 * a;
        std::vector<double> sg(n), raw(n);
        for (int i = 0; i < n; i++) {
            double s = s_max * i / (n - 1.0);
            sg[i] = s;
            int M = 1000; double h = zmax / M, sum = 0.0;
            for (int k = 0; k <= M; k++) {
                double zp = k * h;
                double r  = std::sqrt(zp * zp + s * s);
                double rho = 1.0 / (1.0 + std::exp((r - RA) / a));
                double w = (k == 0 || k == M) ? 1.0 : (k % 2 ? 4.0 : 2.0);
                sum += w * rho;
            }
            raw[i] = 2.0 * sum * h / 3.0;   // 2x: line-of-sight integral is symmetric in z'
        }
        double h2 = sg[1] - sg[0], sum2 = 0.0;
        for (int i = 0; i < n; i++) {
            double w = (i == 0 || i == n - 1) ? 1.0 : (i % 2 ? 4.0 : 2.0);
            sum2 += w * raw[i] * sg[i];
        }
        double integral = 2.0 * M_PI * sum2 * h2 / 3.0;
        std::vector<double> Tg(n);
        for (int i = 0; i < n; i++) Tg[i] = raw[i] * (B_mass / integral);

        spline = gsl_spline_alloc(gsl_interp_cspline, n);
        acc = gsl_interp_accel_alloc();
        gsl_spline_init(spline, sg.data(), Tg.data(), n);
    }
    WSThickness(const WSThickness&) = delete;
    WSThickness& operator=(const WSThickness&) = delete;
    double s_max() const { return s_max_; }
    double operator()(double s) const { return (s >= s_max_) ? 0.0 : gsl_spline_eval(spline, s, acc); }
    ~WSThickness() { gsl_spline_free(spline); gsl_interp_accel_free(acc); }
};

// H_channel(r) = \int d^2s T_B(s) Gamma_AB(|r-s|) P_EMD^channel(|r-s|)
// = \int_0^{s_max} s ds T_B(s) [\int_0^{2pi} dphi Gamma_AB(b(r,s,phi)) P_EMD(b)]
// No extra Theta(b >= bmin) cut: Gamma_AB already removes the overlapping
// configurations (as in Eq. 4-6 of arXiv:2404.09731).
// Built once as an r-spline over [0, r_switch] (each grid point needs the
// (s,phi) double integral, the expensive part); for r >= r_switch (i.e. once
// r is many nuclear thicknesses away, s becomes negligible next to r), the
// exact asymptotic form B_mass*Gamma_AB(r)*P_EMD(r) is used directly instead
// of extending the grid -- continuity across r_switch is a good self-check
// (see scan_flux_eff.cpp, where this was validated before being productionized
// here).
class EffFluxRadial {
private:
    gsl_spline* spline;
    gsl_interp_accel* acc;
    double r_min_, r_switch_;
    double B_mass, S;
    std::string channel;

    static double channel_factor(const std::string& channel, double b, double S)
    {
        double P_no_emd = std::exp(-S / (b * b));
        if (channel == "An0n") return P_no_emd;
        if (channel == "Xn0n") return P_no_emd * (1.0 - P_no_emd);
        if (channel == "0n0n") return P_no_emd * P_no_emd;
        return 1.0;   // AnAn
    }

    static double H_raw(double r, const WSThickness& TB, const std::string& channel, double S)
    {
        int Ns = 150, Nphi = 100;
        double sh = TB.s_max() / Ns, sum_s = 0.0;
        for (int i = 0; i <= Ns; i++) {
            double s = i * sh;
            double Tval = TB(s);
            double ph = M_PI / Nphi, sum_phi = 0.0;
            for (int j = 0; j <= Nphi; j++) {
                double phi = j * ph;
                double b = std::sqrt(std::max(0.0, r * r + s * s - 2.0 * r * s * std::cos(phi)));   // max: rounding at r=s, phi=0
                double val = gamma_aa()(b) * channel_factor(channel, b, S);
                double w = (j == 0 || j == Nphi) ? 1.0 : (j % 2 ? 4.0 : 2.0);
                sum_phi += w * val;
            }
            double phi_integral = 2.0 * (sum_phi * ph / 3.0);   // x2: [0,2pi] via symmetry about phi=pi
            double w = (i == 0 || i == Ns) ? 1.0 : (i % 2 ? 4.0 : 2.0);
            sum_s += w * s * Tval * phi_integral;
        }
        return sum_s * sh / 3.0;
    }

public:
    // r-grid is log-spaced from r0 to r_switch: H(r) varies fastest near
    // R_A (~33 GeV^-1) and needs resolution there, not just far out where
    // it's flat. r_switch itself was found to need pushing well past R_A
    // for a clean match to the asymptotic formula: the EMD-bearing channels
    // (An0n, Xn0n) still carry a percent-level s/r correction at r=300
    // (s reaches up to s_max~66, not yet negligible against r there), which
    // showed up as a visible bump/dip in f_eff(z) right where a given z's
    // r_max(z) integration bound crosses r_switch (see photon_flux_discrepancy.pdf).
    EffFluxRadial(const WSThickness& TB, const std::string& channel_, double S_,
                  double B_mass_, double r_switch, int nr = 400, double r0 = 0.5)
        : r_min_(r0), r_switch_(r_switch), B_mass(B_mass_), S(S_), channel(channel_)
    {
        std::vector<double> lrg(nr), Hg(nr);
        double lr0 = std::log(r0), lr1 = std::log(r_switch);
        for (int i = 0; i < nr; i++) {
            double lr = lr0 + (lr1 - lr0) * i / (nr - 1.0);
            lrg[i] = lr;
            Hg[i] = H_raw(std::exp(lr), TB, channel, S);
        }
        spline = gsl_spline_alloc(gsl_interp_cspline, nr);
        acc = gsl_interp_accel_alloc();
        gsl_spline_init(spline, lrg.data(), Hg.data(), nr);
    }
    EffFluxRadial(const EffFluxRadial&) = delete;
    EffFluxRadial& operator=(const EffFluxRadial&) = delete;

    double operator()(double r) const
    {
        if (r >= r_switch_) return B_mass * gamma_aa()(r) * channel_factor(channel, r, S);
        double rc = (r < r_min_) ? r_min_ : r;
        return gsl_spline_eval(spline, std::log(rc), acc);
    }

    ~EffFluxRadial() { gsl_spline_free(spline); gsl_interp_accel_free(acc); }
};


void load_data_and_initialize(const std::string& filename);

// Builds the Woods-Saxon form factor spline (module-level, like gamma_aa()
// below). Must be called before flux_density_WS()/photon_flux() with
// par->flux_model=="WS" is used. Safe to call more than once (rebuilds).
void init_ws_form_factor(double RA_fm = 6.49, double a_fm = 0.54, double mn = 0.938);

// Builds the "effective flux" table (Eq. 4 of arXiv:2404.09731) for the
// given channel, replacing the b-integral entirely: effective_photon_flux()
// below then depends only on qp (equivalently z), not on b at all. Requires
// init_ws_form_factor() to already have been called (the bare flux used
// inside the r-integral is the WS one, needed since r ranges down to 0 here
// -- see EffFluxRadial). par is only used for alpha/Z/mn/ss/S (same
// struct as flux_density_WS()).
void init_effective_flux(const std::string& channel, void* par,
                          double RA_fm = 6.49, double a_fm = 0.54, double B_mass = 208.0);

// FLUX_MODEL=STARLIGHT: P. Paakkinen's tabulated WS effective flux
// (inputs/photon_flux/log-flux-tbl-WS.dta, arXiv:2404.09731; sqrt(s_NN) =
// 5.36 TeV, sigma_NN = 90.8533 mb, R_WS = 6.49 fm, d_WS = 0.54 fm), read
// once here. Only AnAn and An0n are tabulated. effective_photon_flux()
// then returns it instead of our own f_eff.
void init_starlight_flux(const std::string& channel,
                         const std::string& filename = "./inputs/photon_flux/log-flux-tbl-WS.dta");

// Effective photon flux, i.e. f_eff(z) with the full b (and target-nucleus
// r,s geometric convolution) already integrated out -- returns dN/domega,
// the same convention as \int db photon_flux(b,qp,par) would give under the
// old single-b treatment, so it's a drop-in replacement for that integral,
// not for photon_flux() itself (there is no b left to pass in).
double effective_photon_flux(double qp, void* par);

GammaAA& gamma_aa();

// Proton-target (pA) hadronic survival factor, Eq. 26 of arXiv:2606.05469:
// Gamma_pA(b) = exp(-sigma_NN * T_A(b)), the optical-limit probability that
// the proton does NOT interact hadronically with the (lead) nucleus emitting
// the photon, at impact parameter b. Unlike Gamma_AA (nucleus-nucleus, read
// from a precomputed Glauber table), this is a simple analytic exponential
// of the SAME Woods-Saxon thickness function T_A used by the effective-flux
// convolution above -- built directly here rather than tabulated, since the
// proton itself contributes no extended structure to convolve over (it's the
// pointlike target in this survival factor, unlike the emitting nucleus).
// par->target=="pA" selects this inside photon_flux() (photon_flux.cpp);
// no EMD suppression is applied for pA (see Sec. 6.1 of that paper).
void init_pA_flux(double sigma_NN_mb, double RA_fm = 6.49, double a_fm = 0.54, double B_mass = 208.0);
double GammaPA(double b);

#endif
