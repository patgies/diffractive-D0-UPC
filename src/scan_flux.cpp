// Diagnostic: b-integrated photon flux from photon_flux() (gamma_aa.cpp),
// as dN/dy with y = 2*omega/sqrt(s_NN), for comparison with the tabulated
// effective fluxes of arXiv:2404.09731 (inputs/photon_flux/, see
// python/flux_comparison.py). Not part of the physics pipeline.
//
// usage: ./build/bin/scan_flux > out/flux_scan.csv
// columns: channel,bmax_GeV-1,y,dN_dy   (dN_dy = dN/domega * sqrt(s)/2)
#include "def.hpp"
#include "gamma_aa.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

int main()
{
    load_data_and_initialize("./inputs/Gamma_AA.dat");

    parameters par;
    par.ss    = 5360.0;
    par.alpha = 1.0/137.0;
    par.Z     = 82.0;
    par.mn    = 0.938;
    par.S     = std::pow(17.4, 2) / std::pow(0.197327, 2);
    const char* bm = std::getenv("BMIN_FM");   // default: same lower cut as main.cpp (14.2 fm)
    par.bmin  = (bm ? std::atof(bm) : 14.2) / 0.197327;

    const std::vector<std::string> channels = {"PL(AnAn)", "AnAn", "An0n"};
    // 650 GeV^-1: the upper b limit actually used in main.cpp; 1e7: effectively unbounded
    const std::vector<double> bmaxs = {650.0, 1.0e7};

    std::cout << "channel,bmax,y,dN_dy\n" << std::setprecision(10);
    const int ny = 60, nb = 6000;
    for (const auto& ch : channels) {
        par.channel = ch;
        for (double bmax : bmaxs) {
            for (int iy = 0; iy < ny; iy++) {
                double y  = std::pow(10.0, -4.0 + 4.0 * iy / (ny - 1.0));
                double qp = y * par.ss / std::sqrt(2.0);      // omega = qp/sqrt(2) = y*sqrt(s)/2
                // integrate photon_flux(b) db on a log-b grid (Simpson in ln b)
                double a = std::log(par.bmin), c = std::log(bmax), h = (c - a) / nb, sum = 0.0;
                for (int i = 0; i <= nb; i++) {
                    double b = std::exp(a + i * h);
                    double w = (i == 0 || i == nb) ? 1.0 : (i % 2 ? 4.0 : 2.0);
                    sum += w * photon_flux(b, qp, &par) * b;   // db = b dlnb
                }
                double dN_domega = sum * h / 3.0;
                double dN_dy = dN_domega * par.ss / 2.0;
                std::cout << ch << "," << bmax << "," << y << "," << dN_dy << "\n";
            }
        }
    }
    return 0;
}
