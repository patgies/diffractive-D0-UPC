// Diagnostic: b-integrated photon flux from photon_flux() (photon_flux.cpp),
// as dN/dy with y = 2*omega/sqrt(s_NN), for comparison with the tabulated
// effective fluxes of arXiv:2404.09731 (inputs/photon_flux/, see
// python/flux_comparison.py). Not part of the physics pipeline.
//
// usage: ./build/bin/scan_flux > out/flux_scan.csv
// columns: channel,y,dN_dy   (dN_dy = dN/domega * sqrt(s)/2)
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
    load_data_and_initialize("./inputs/Gamma_AA.dat");

    parameters par;
    par.ss    = 5360.0;   // GeV, sqrt(s_NN) of the Starlight tables (inputs/photon_flux/*.dta headers)
    par.alpha = 1.0/137.0;
    par.Z     = 82.0;
    par.mn    = 0.938;
    par.S     = std::pow(17.4, 2) / std::pow(0.197327, 2);
    // default: b_min = 14.205 fm, matched to the arXiv:2404.09731 PL table
    // (ratio = 1 to <0.1% up to y = 0.6) -- 14.2 fm (main.cpp) gives an
    // exp(0.047 y) excess over the table at large y, 72 GeV^-1 a deficit.
    const char* bm = std::getenv("BMIN_FM");
    par.bmin  = (bm ? std::atof(bm) : 14.205) / 0.197327;
    par.flux_model = std::getenv("FLUX_MODEL") ? std::getenv("FLUX_MODEL") : "PL";
    if (par.flux_model == "WS") init_ws_form_factor();

    const std::vector<std::string> channels = {"PL(AnAn)", "AnAn", "An0n", "Xn0n"};

    std::cout << "channel,y,dN_dy\n" << std::setprecision(10);
    const int ny = 60, nb = 6000;
    for (const auto& ch : channels) {
        par.channel = ch;
        for (int iy = 0; iy < ny; iy++) {
            double y    = std::pow(10.0, -4.0 + 4.0 * iy / (ny - 1.0));
            double qp   = y * par.ss / std::sqrt(2.0);    // omega = qp/sqrt(2) = y*sqrt(s)/2
            // bmax = 60/(y*mn): large enough that eta=y*mn*b comfortably
            // clears the eta>50 cutoff already built into photon_flux()
            // (photon_flux.cpp), so this reproduces the true b->infinity
            // integral at every y without wasting effort far past where the
            // flux is already zero. Same choice now used for par.bmax in
            // main.cpp/main_xpom.cpp (see there for why it has to be
            // computed differently there: q+/y is a Monte Carlo variable in
            // the main pipeline, not looped over externally like here).
            double bmax = 60.0 / (y * par.mn);
            // integrate photon_flux(b) db on a log-b grid (Simpson in ln b)
            // sharp bmin cut only for the point-like channel; the others
            // start (effectively) at b = 0, Gamma_AA does the suppression
            double b_lo = (ch == "PL(AnAn)") ? par.bmin : 1e-3;
            double a = std::log(b_lo), c = std::log(bmax), h = (c - a) / nb, sum = 0.0;
            for (int i = 0; i <= nb; i++) {
                double b = std::exp(a + i * h);
                double w = (i == 0 || i == nb) ? 1.0 : (i % 2 ? 4.0 : 2.0);
                sum += w * photon_flux(b, qp, &par) * b;   // db = b dlnb
            }
            double dN_domega = sum * h / 3.0;
            double dN_dy = dN_domega * par.ss / 2.0;
            std::cout << ch << "," << y << "," << dN_dy << "\n";
        }
    }
    return 0;
}
