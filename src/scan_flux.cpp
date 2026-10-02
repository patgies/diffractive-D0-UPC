// Test tool: photon flux dN/dy integrated over b, for python/flux_comparison.py.
// usage: ./build/bin/scan_flux > output/flux_scan/flux_scan.dat
#include "def.hpp"
#include "photon_flux.hpp"
#include <gsl/gsl_errno.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

int main()
{
    gsl_set_error_handler_off();
    load_data_and_initialize("./inputs/WS_photon_flux/Gamma_AA.dat");

    parameters par;
    par.ss    = 5360.0;   // GeV, sqrt(s_NN) of the Starlight tables (inputs/Starlight_photon_flux/*.dta headers)
    par.alpha = 1.0/137.0;
    par.Z     = 82.0;
    par.mn    = 0.938;
    par.S     = std::pow(17.4, 2) / std::pow(0.197327, 2);
    // default b_min = 14.205 fm, the value that agrees with the PL table of arXiv:2404.09731
    const char* bm = std::getenv("BMIN_FM");
    par.bmin  = (bm ? std::atof(bm) : 14.205) / 0.197327;
    par.flux_model = std::getenv("FLUX_MODEL") ? std::getenv("FLUX_MODEL") : "PL";
    if (par.flux_model == "WS") init_ws_form_factor();

    const std::vector<std::string> channels = {"PL(AnAn)", "AnAn", "An0n", "Xn0n"};

    std::cout << "# channel  y  dN_dy\n" << std::setprecision(10);
    const int ny = 60, nb = 6000;
    for (const auto& ch : channels) {
        par.channel = ch;
        for (int iy = 0; iy < ny; iy++) {
            double y    = std::pow(10.0, -4.0 + 4.0 * iy / (ny - 1.0));
            double qp   = y * par.ss / std::sqrt(2.0);    // omega = qp/sqrt(2) = y*sqrt(s)/2
            // bmax = 60/(y*mn): beyond eta > 50, where photon_flux() is already zero
            double bmax = 60.0 / (y * par.mn);
            // Simpson rule in ln b. Only the point-like channel starts at bmin
            double b_lo = (ch == "PL(AnAn)") ? par.bmin : 1e-3;
            double a = std::log(b_lo), c = std::log(bmax), h = (c - a) / nb, sum = 0.0;
            for (int i = 0; i <= nb; i++) {
                double b = std::exp(a + i * h);
                double w = (i == 0 || i == nb) ? 1.0 : (i % 2 ? 4.0 : 2.0);
                sum += w * photon_flux(b, qp, &par) * b;   // db = b dlnb
            }
            double dN_domega = sum * h / 3.0;
            double dN_dy = dN_domega * par.ss / 2.0;
            std::cout << ch << "  " << y << "  " << dN_dy << "\n";
        }
    }
    return 0;
}
