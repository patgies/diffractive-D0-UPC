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


// Woods-Saxon charge form factor F_A(q), formula of Maximon & Schrack,
// J. Res. Natl. Bur. Stand. B 70, 85 (1966). F(0) = 1.
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
    // RA and a in GeV^-1, q_max in GeV.
    WSFormFactor(double RA, double a, double q_max) : c_(RA), a_(a), q_max_(q_max)
    {
        double bracket = (4.0 * M_PI * c_ / 3.0) * (M_PI * a_ * M_PI * a_ + c_ * c_)
                        + 8.0 * M_PI * a_ * a_ * a_ * series_cube(c_, a_);
        rho0_ = 1.0 / bracket;
    }

    WSFormFactor(const WSFormFactor&) = delete;
    WSFormFactor& operator=(const WSFormFactor&) = delete;

    double operator()(double q) const
    {
        if (q >= q_max_) return 0.0;      // F is already zero before q_max
        if (q < 1e-6) return 1.0;         // F(0)=1. This avoids 0/0 in the formula below
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

// Table of the k_perp integral I(z,b) of the Woods-Saxon flux, on a grid in (log z, log b).
// f_WS(z,b) = alpha Z^2/(pi^2 z) * I(z,b)^2 (arXiv:2404.09731 Eq. 11).
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
        std::vector<double> grid_z(nz * nb);   // 1D array in the order gsl_spline2d needs
        for (int i = 0; i < nz; i++) {
            double z = std::exp(lz[i]);
            for (int j = 0; j < nb; j++) {
                double b = std::exp(lb[j]);
                double I = raw_I(z, b, mn, FA, q_max);
                // if the k_perp integral fails (it can near z=1), store a tiny value so the spline still works
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

    // Returns I(z,b)^2, without the alpha Z^2/(pi^2 z) prefactor.
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
 
        // keep the value between 0 and 1: the spline can give ~-1e-32 where Gamma -> 0
        return std::min(1.0, std::max(0.0, TA(b)));
    }
};

GammaAA& gamma_aa();   // declared here so EffFluxRadial below can use it

// Nuclear thickness function T_B(s) from the Woods-Saxon density,
// scaled so that it adds up to B_mass nucleons.
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
            raw[i] = 2.0 * sum * h / 3.0;   // times 2: the integral is the same for z' > 0 and z' < 0
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

// H_channel(r) = \int d^2s T_B(s) Gamma_AB(|r-s|) P_EMD(|r-s|), stored as a spline in r up to
// r_switch. For larger r we use B_mass*Gamma_AB(r)*P_EMD(r).
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
                double b = std::sqrt(std::max(0.0, r * r + s * s - 2.0 * r * s * std::cos(phi)));   // max avoids a tiny negative number at r=s, phi=0
                double val = gamma_aa()(b) * channel_factor(channel, b, S);
                double w = (j == 0 || j == Nphi) ? 1.0 : (j % 2 ? 4.0 : 2.0);
                sum_phi += w * val;
            }
            double phi_integral = 2.0 * (sum_phi * ph / 3.0);   // times 2: phi from pi to 2pi gives the same as 0 to pi
            double w = (i == 0 || i == Ns) ? 1.0 : (i % 2 ? 4.0 : 2.0);
            sum_s += w * s * Tval * phi_integral;
        }
        return sum_s * sh / 3.0;
    }

public:
    // The r grid is evenly spaced in log r. r_switch must be much larger than R_A (see EFF_R_SWITCH).
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

// Builds the Woods-Saxon flux table. Call it before using flux_model "WS".
void init_ws_form_factor(double RA_fm = 6.49, double a_fm = 0.54, double mn = 0.938);

// Builds the effective flux f_eff(z) for the given channel (arXiv:2404.09731 Eq. 4).
void init_effective_flux(const std::string& channel, void* par,
                          double RA_fm = 6.49, double a_fm = 0.54, double B_mass = 208.0);

// FLUX_MODEL=STARLIGHT: reads the effective flux table of P. Paakkinen
// (arXiv:2404.09731). The table only has AnAn and An0n.
void init_starlight_flux(const std::string& channel,
                         const std::string& filename = "./inputs/Starlight_photon_flux/log-flux-tbl-WS.dta");

// Effective flux dN/domega. The sum over b is already done.
double effective_photon_flux(double qp, void* par);

GammaAA& gamma_aa();

// Proton target: Gamma_pA(b) = exp(-sigma_NN * T_A(b)),
// Eq. 26 of arXiv:2606.05469. There is no EMD factor for pA.
void init_pA_flux(double sigma_NN_mb, double RA_fm = 6.49, double a_fm = 0.54, double B_mass = 208.0);
double GammaPA(double b);

#endif
