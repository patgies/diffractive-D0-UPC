// Diagnostic ONLY (not part of the physics pipeline): tests whether the
// "effective flux" geometric convolution --
//   f_eff(z) = (1/B) \int d2r \int d2s  f_gamma/A(z,r) T_B(s) Gamma_AB(r-s)
// (Eq. 4 of arXiv:2404.09731, Eskola et al.) -- closes the gap between this
// code's simple single-b treatment (par.bmin cut + Gamma_AA(b) at the SAME b
// as the photon flux, i.e. Eq. 1/15 of our own paper and of Guzey et al.
// 2606.05469) and Paakkinen's tabulated "Starlight" flux
// (inputs/photon_flux/log-flux-tbl-WS.dta).
//
// Structure (see the derivation notes inline): building this bottom-up,
//   T_B(s)          -- nuclear thickness function, normalized to integrate
//                       to B=208 nucleons (WSThickness below)
//   G_channel(r,s)  = \int_0^{2pi} dphi Gamma_AB(|r-s|) P_EMD^channel(|r-s|)
//                       (fixed-grid Simpson in phi; no extra b >= bmin cut --
//                       Gamma_AB already removes the overlapping configurations)
//   H_channel(r)    = \int s ds T_B(s) G_channel(r,s)  (fixed-grid Simpson
//                       in s), which asymptotes to B*Gamma_AB(r)*P_EMD(r)
//                       once r is well past where T_B(s) has support --
//                       used as a fast/exact fallback for large r instead of
//                       extending the (s,phi) grid out there.
//   f_eff(z)        = (2*pi/B) \int r dr f_gamma/A^WS(z,r) H_channel(r)
//                       (fixed-grid Simpson in r, log-spaced), using the
//                       already-validated WS bare flux for f_gamma/A(z,r) --
//                       needed here (unlike the main pipeline) because r is
//                       NOT bounded away from 0: the photon can be emitted
//                       from anywhere in nucleus A, including deep inside it.
//
// usage: ./build/bin/scan_flux_eff <channel> > out/flux_scan_eff_<channel>.csv
#include "def.hpp"
#include "photon_flux.hpp"
#include <gsl/gsl_errno.h>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

static const double hbarc = 0.197327;

// channel-dependent electromagnetic-dissociation factor, same formula as
// photon_flux() in photon_flux.cpp
static double channel_emd_factor(const string& channel, double b, double S)
{
    double P_emd    = S / (b * b);
    double P_no_emd = std::exp(-P_emd);
    if (channel == "An0n") return P_no_emd;
    if (channel == "Xn0n") return P_no_emd * (1.0 - P_no_emd);
    return 1.0;   // AnAn
}

// T_B(s): nuclear thickness function (line-of-sight integral of the
// Woods-Saxon density), normalized so that 2*pi*\int T_B(s) s ds = B.
struct Thickness {
    gsl_spline* spline;
    gsl_interp_accel* acc;
    double s_max;

    Thickness(double RA, double a, double B_mass, double s_max_, int n = 150) : s_max(s_max_)
    {
        double zmax = RA + 12.0 * a;
        std::vector<double> sg(n), raw(n);
        for (int i = 0; i < n; i++) {
            double s = s_max * i / (n - 1.0);
            sg[i] = s;
            // 2*\int_0^zmax dz' rho(sqrt(z'^2+s^2)), fixed Simpson (smooth integrand)
            int M = 2000; double h = zmax / M, sum = 0.0;
            for (int k = 0; k <= M; k++) {
                double zp = k * h;
                double r  = std::sqrt(zp * zp + s * s);
                double rho = 1.0 / (1.0 + std::exp((r - RA) / a));
                double w = (k == 0 || k == M) ? 1.0 : (k % 2 ? 4.0 : 2.0);
                sum += w * rho;
            }
            raw[i] = 2.0 * sum * h / 3.0;
        }
        // normalize: 2*pi*\int raw(s) s ds == B_mass
        int n2 = n; double h2 = sg[1] - sg[0], sum2 = 0.0;
        for (int i = 0; i < n2; i++) {
            double w = (i == 0 || i == n2 - 1) ? 1.0 : (i % 2 ? 4.0 : 2.0);
            sum2 += w * raw[i] * sg[i];
        }
        double integral = 2.0 * M_PI * sum2 * h2 / 3.0;
        std::vector<double> Tg(n);
        for (int i = 0; i < n; i++) Tg[i] = raw[i] * (B_mass / integral);

        spline = gsl_spline_alloc(gsl_interp_cspline, n);
        acc = gsl_interp_accel_alloc();
        gsl_spline_init(spline, sg.data(), Tg.data(), n);
    }
    double operator()(double s) const { return (s >= s_max) ? 0.0 : gsl_spline_eval(spline, s, acc); }
    ~Thickness() { gsl_spline_free(spline); gsl_interp_accel_free(acc); }
};

// H_channel(r) = \int_0^{s_max} s ds T_B(s) * [\int_0^{2pi} dphi Gamma_AB(|r-s|) P_EMD(|r-s|)]
// Exact asymptotic fallback for r >= r_switch: H(r) -> B*Gamma_AB(r)*P_EMD(r).
static double H_channel(double r, const Thickness& TB, const string& channel, double S,
                         double B_mass, double r_switch)
{
    if (r >= r_switch) {
        return B_mass * gamma_aa()(r) * channel_emd_factor(channel, r, S);
    }
    // grid sizes, overridable for convergence tests (NS, NPHI env vars)
    static const int Ns   = std::getenv("NS")   ? std::atoi(std::getenv("NS"))   : 120;
    static const int Nphi = std::getenv("NPHI") ? std::atoi(std::getenv("NPHI")) : 80;
    double sh = TB.s_max / Ns, sum_s = 0.0;
    for (int i = 0; i <= Ns; i++) {
        double s = i * sh;
        double Tval = TB(s);
        if (Tval <= 0.0 && i > 0) continue;
        double ph = M_PI / Nphi, sum_phi = 0.0;
        for (int j = 0; j <= Nphi; j++) {
            double phi = j * ph;
            double b = std::sqrt(std::max(0.0, r * r + s * s - 2.0 * r * s * std::cos(phi)));   // max: rounding at r=s, phi=0
            double val = gamma_aa()(b) * channel_emd_factor(channel, b, S);
            double w = (j == 0 || j == Nphi) ? 1.0 : (j % 2 ? 4.0 : 2.0);
            sum_phi += w * val;
        }
        double phi_integral = 2.0 * (sum_phi * ph / 3.0);   // *2 for [0,2pi] via symmetry about phi=pi
        double w = (i == 0 || i == Ns) ? 1.0 : (i % 2 ? 4.0 : 2.0);
        sum_s += w * s * Tval * phi_integral;
    }
    return sum_s * sh / 3.0;
}

int main(int argc, char** argv)
{
    gsl_set_error_handler_off();
    string channel = (argc > 1) ? argv[1] : "AnAn";

    load_data_and_initialize("./inputs/Gamma_AA.dat");
    init_ws_form_factor();   // builds the WS bare-flux table (needed for r < R_A)

    parameters par;
    par.ss    = 5360.0;   // GeV, sqrt(s_NN) of the Starlight tables (inputs/photon_flux/*.dta headers)
    par.alpha = 1.0 / 137.0;
    par.Z     = 82.0;
    par.mn    = 0.938;
    par.S     = std::pow(17.4, 2) / std::pow(0.197327, 2);
    par.bmin  = 14.2 / 0.197327;
    par.flux_model = "WS";

    double RA = 6.49 / hbarc, a = 0.54 / hbarc;
    double B_mass = 208.0;
    double s_max = RA + 12.0 * a;
    double r_switch = 800.0;   // GeV^-1; well past where T_B(s) has support -- 300 left a
                               // percent-level s/r mismatch for the EMD-bearing channels

    Thickness TB(RA, a, B_mass, s_max);

    // Sanity check: as r -> large, H_channel should match the asymptotic
    // formula smoothly (continuity check across r_switch).
    {
        double r1 = r_switch * 0.98, r2 = r_switch * 1.02;
        double h1 = H_channel(r1, TB, channel, par.S, B_mass, r_switch);
        double h2 = H_channel(r2, TB, channel, par.S, B_mass, r_switch);
        std::cerr << "# continuity check at r_switch=" << r_switch
                  << ": H(" << r1 << ")=" << h1 << "  H(" << r2 << ")=" << h2
                  << "  ratio=" << h1/h2 << "\n";
    }

    std::cout << "channel,y,dN_dy\n" << std::setprecision(10);
    const int ny = 120;
    for (int iy = 0; iy < ny; iy++) {
        double y  = std::pow(10.0, -4.0 + 4.0 * iy / (ny - 1.0));   // z = photon energy fraction
        double qp = y * par.ss / std::sqrt(2.0);
        double r_max = 60.0 / (y * par.mn);   // same adaptive cutoff idea as scan_flux.cpp

        // r-integral: log-spaced from a small r0 up to r_max, fixed Simpson in ln r
        double r0 = 1e-3;
        static const int Nr = std::getenv("NR") ? std::atoi(std::getenv("NR")) : 4000;
        double a_ = std::log(r0), c_ = std::log(r_max), h_ = (c_ - a_) / Nr, sum = 0.0;
        for (int i = 0; i <= Nr; i++) {
            double r = std::exp(a_ + i * h_);
            double flux = flux_density_WS(y, r, &par);   // bare flux, valid down to r=0
            double Hval = H_channel(r, TB, channel, par.S, B_mass, r_switch);
            double w = (i == 0 || i == Nr) ? 1.0 : (i % 2 ? 4.0 : 2.0);
            sum += w * (r * flux * Hval) * r;   // extra r: dr = r dlnr
        }
        double r_integral = sum * h_ / 3.0;
        double f_eff = (2.0 * M_PI / B_mass) * r_integral;

        double dN_domega = f_eff;               // f_eff(z) already IS dN/domega-like (matches flux_density_WS's own convention)
        double dN_dy = dN_domega * par.ss / 2.0;
        std::cout << channel << "," << y << "," << dN_dy << "\n";
    }
    return 0;
}
