// Test tool: effective flux f_eff(z) (Eq. 4 of arXiv:2404.09731) for one channel.
// usage: ./build/bin/scan_flux_eff <channel> > output/flux_scan/flux_scan_eff_<channel>.dat
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

// EMD factor of each channel, same formula as in photon_flux().
static double channel_emd_factor(const string& channel, double b, double S)
{
    double P_emd    = S / (b * b);
    double P_no_emd = std::exp(-P_emd);
    if (channel == "An0n") return P_no_emd;
    if (channel == "Xn0n") return P_no_emd * (1.0 - P_no_emd);
    if (channel == "0n0n") return P_no_emd * P_no_emd;
    return 1.0;   // AnAn
}

// T_B(s): nuclear thickness function, scaled so that it adds up to B nucleons.
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
            // integral of the density along z (Simpson rule)
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

// H_channel(r) = \int s ds T_B(s) \int dphi Gamma_AB(|r-s|) P_EMD(|r-s|);
// for r >= r_switch we use B*Gamma_AB(r)*P_EMD(r).
static double H_channel(double r, const Thickness& TB, const string& channel, double S,
                         double B_mass, double r_switch)
{
    if (r >= r_switch) {
        return B_mass * gamma_aa()(r) * channel_emd_factor(channel, r, S);
    }
    // grid sizes, can be changed for tests with NS and NPHI
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
            double b = std::sqrt(std::max(0.0, r * r + s * s - 2.0 * r * s * std::cos(phi)));   // max avoids a tiny negative number at r=s, phi=0
            double val = gamma_aa()(b) * channel_emd_factor(channel, b, S);
            double w = (j == 0 || j == Nphi) ? 1.0 : (j % 2 ? 4.0 : 2.0);
            sum_phi += w * val;
        }
        double phi_integral = 2.0 * (sum_phi * ph / 3.0);   // times 2: phi from pi to 2pi gives the same as 0 to pi
        double w = (i == 0 || i == Ns) ? 1.0 : (i % 2 ? 4.0 : 2.0);
        sum_s += w * s * Tval * phi_integral;
    }
    return sum_s * sh / 3.0;
}

int main(int argc, char** argv)
{
    gsl_set_error_handler_off();
    string channel = (argc > 1) ? argv[1] : "AnAn";

    load_data_and_initialize("./inputs/WS_photon_flux/Gamma_AA.dat");
    init_ws_form_factor();   // builds the WS flux table (needed for r < R_A)

    parameters par;
    par.ss    = 5360.0;   // GeV, sqrt(s_NN) of the Starlight tables (inputs/Starlight_photon_flux/*.dta headers)
    par.alpha = 1.0 / 137.0;
    par.Z     = 82.0;
    par.mn    = 0.938;
    par.S     = std::pow(17.4, 2) / std::pow(0.197327, 2);
    par.bmin  = 14.2 / 0.197327;
    par.flux_model = "WS";

    double RA = 6.49 / hbarc, a = 0.54 / hbarc;
    double B_mass = 208.0;
    double s_max = RA + 12.0 * a;
    double r_switch = 800.0;   // GeV^-1, much larger than the size of the nucleus

    Thickness TB(RA, a, B_mass, s_max);

    // check that H_channel does not jump at r_switch
    {
        double r1 = r_switch * 0.98, r2 = r_switch * 1.02;
        double h1 = H_channel(r1, TB, channel, par.S, B_mass, r_switch);
        double h2 = H_channel(r2, TB, channel, par.S, B_mass, r_switch);
        std::cerr << "# continuity check at r_switch=" << r_switch
                  << ": H(" << r1 << ")=" << h1 << "  H(" << r2 << ")=" << h2
                  << "  ratio=" << h1/h2 << "\n";
    }

    std::cout << "# channel  y  dN_dy\n" << std::setprecision(10);
    const int ny = 120;
    for (int iy = 0; iy < ny; iy++) {
        double y  = std::pow(10.0, -4.0 + 4.0 * iy / (ny - 1.0));   // z = photon energy fraction
        double r_max = 60.0 / (y * par.mn);   // same upper limit as in scan_flux.cpp

        // r integral from r0 to r_max, Simpson rule in ln r
        double r0 = 1e-3;
        static const int Nr = std::getenv("NR") ? std::atoi(std::getenv("NR")) : 4000;
        double a_ = std::log(r0), c_ = std::log(r_max), h_ = (c_ - a_) / Nr, sum = 0.0;
        for (int i = 0; i <= Nr; i++) {
            double r = std::exp(a_ + i * h_);
            double flux = flux_density_WS(y, r, &par);   // WS flux, works down to r=0
            double Hval = H_channel(r, TB, channel, par.S, B_mass, r_switch);
            double w = (i == 0 || i == Nr) ? 1.0 : (i % 2 ? 4.0 : 2.0);
            sum += w * (r * flux * Hval) * r;   // extra r: dr = r dlnr
        }
        double r_integral = sum * h_ / 3.0;
        double f_eff = (2.0 * M_PI / B_mass) * r_integral;

        double dN_domega = f_eff;               // f_eff(z) is already dN/domega (same units as flux_density_WS)
        double dN_dy = dN_domega * par.ss / 2.0;
        std::cout << channel << "  " << y << "  " << dN_dy << "\n";
    }
    return 0;
}
