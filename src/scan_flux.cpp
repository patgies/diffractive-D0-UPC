#include "def.hpp"
#include "photon_flux.hpp"
#include <gsl/gsl_errno.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

// dN/domega from the flux at one b, integrated over b (Simpson rule in ln b)
static double flux_over_b(double z_gamma, double qp, parameters& par)
{
    const int nb = 6000;
    double bmax = 60.0 / (z_gamma * par.mn);   // beyond eta > 50 photon_flux() is already zero
    double b_lo = (par.channel == "PL(AnAn)" || par.target == "pA") ? par.bmin : 1e-3;   // see b_integration_min()
    double a = std::log(b_lo), c = std::log(bmax), h = (c - a) / nb, sum = 0.0;
    for (int i = 0; i <= nb; i++) {
        double b = std::exp(a + i * h);
        double w = (i == 0 || i == nb) ? 1.0 : (i % 2 ? 4.0 : 2.0);
        sum += w * photon_flux(b, qp, &par) * b;   // db = b dlnb
    }
    return sum * h / 3.0;
}

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <EFF|STARLIGHT|PL|WS> <channel>...   (TARGET=pA: WS or PL, channel pA)" << std::endl;
        return 1;
    }
    gsl_set_error_handler_off();
    const std::string model = argv[1];

    load_data_and_initialize("./input/WS_photon_flux/Gamma_AA.dat");

    parameters par;
    par.ss    = 5360.0;   // GeV, sqrt(s_NN) of the Starlight tables (input/Starlight_photon_flux/*.dta headers)
    par.alpha = 1.0 / 137.0;
    par.Z     = 82.0;
    par.mn    = 0.938;
    par.S     = std::pow(17.4, 2) / std::pow(0.197327, 2);
    // b_min = 14.205 fm, the value that agrees with the PL table of arXiv:2404.09731
    const char* bm = std::getenv("BMIN_FM");
    par.bmin  = (bm ? std::atof(bm) : 14.205) / 0.197327;
    par.flux_model = model;
    if (model == "WS") init_ws_form_factor();

    // TARGET=pA: flux of the lead nucleus in p+Pb at 8.16 TeV, with Gamma_pA(b) and no EMD (as in main.cpp)
    const char* target = std::getenv("TARGET");
    if (target && std::string(target) == "pA") {
        if (model != "WS" && model != "PL") {
            std::cerr << "Error: TARGET=pA needs WS or PL." << std::endl;
            return 1;
        }
        par.target = "pA";
        par.ss     = 8160.0;
        double sigma_NN_mb = std::getenv("SIGMA_NN") ? std::atof(std::getenv("SIGMA_NN")) : 99.0;
        init_pA_flux(sigma_NN_mb);
        const double R_PL_fm = 7.1;   // point-like radius (see main.cpp)
        par.bmin = (model == "PL") ? ((bm ? std::atof(bm) : 1.1 * R_PL_fm) / 0.197327) : 1e-3;   // BMIN_FM changes it
        std::cout << "# target         : pA, sqrt(s_NN) = " << par.ss << " GeV, sigma_NN = " << sigma_NN_mb << " mb\n";
    }

    std::cout << "# flux_model     : " << model << "\n";
    std::cout << "# Gamma_AA table : " << gamma_aa_info() << "\n";
    std::cout << "# channel  y  dN_dy\n" << std::setprecision(10);
    const int ny = 120;
    for (int k = 2; k < argc; k++) {
        par.channel = argv[k];
        if (model == "EFF")       init_effective_flux(par.channel, &par);
        if (model == "STARLIGHT") init_starlight_flux(par.channel);
        for (int iy = 0; iy < ny; iy++) {
            double z_gamma = std::pow(10.0, -4.0 + 4.0 * iy / (ny - 1.0));   // photon energy fraction
            double qp      = z_gamma * par.ss / std::sqrt(2.0);   // omega = qp/sqrt(2) = z_gamma*sqrt(s)/2
            double dN_domega = is_effective_flux(&par) ? effective_photon_flux(qp, &par) : flux_over_b(z_gamma, qp, par);
            std::cout << par.channel << "  " << z_gamma << "  " << dN_domega * par.ss / 2.0 << "\n";
        }
    }
    return 0;
}
