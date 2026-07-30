#include "amplitudelib.hpp"
#include "def.hpp"
#include "tools.hpp"
#include "gamma_aa.hpp"
#include <string>
#include <iostream>
#include <cmath>
#include <ctime>
#include <chrono>
#include <vector>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>

using namespace std;

// Fixed-W (no photon-flux) exclusive and diffractive spectra, as a function
// of the observed quark's transverse momentum K, at a fixed rapidity y.
// "Fixed W" means we set the photon energy qp to 3*pp by hand (so z1=1/3).

int main(int argc, char* argv[])
{
    string datafile = "./data/proton/mve.dat";

    AmplitudeLib inst(datafile);
    inst.SetOutOfRangeErrors(false);

    gsl_set_error_handler_off();
    gsl_rng_env_setup();

    load_data_and_initialize("./data/Gamma_AA.dat");   // unused (no flux here)


    parameters param;
    param.dipole   = &inst;
    param.datafile = datafile;
    param.ss       = 5360.0;
    param.calls_diff = (argc > 3) ? (size_t)atof(argv[3]) : (size_t)1e5;
    param.calls_excl = param.calls_diff;
    param.m        = 1.5;
    param.m2       = param.m * param.m;
    param.y        = (argc > 1) ? atof(argv[1]) : 0.0;
    param.kmax     = 400.0;

    int n_pt = (argc > 2) ? atoi(argv[2]) : 40;
    double pt_min = 0.2, pt_max = 20.0;

    std::time_t now = std::time(nullptr);
    char date_buf[32];
    std::strftime(date_buf, sizeof(date_buf), "%Y-%m-%d %H:%M:%S", std::localtime(&now));

    cout << "# ============================================================" << endl;
    cout << "# Date        =  " << date_buf << endl;
    cout << "# Fixed-W (no photon flux) exclusive & diffractive dsigma/d2K" << endl;
    cout << "# W^2 = sqrt(2)*ss*3*pp fixed via qp=3*pp (z1=1/3); y fixed" << endl;
    cout << "# x_pom > 0.1 forced to 0 on both sides" << endl;
    cout << "# exclusive column:   missing prefactor alpha*Nc*ef^2*Sperp/(2pi^2)" << endl;
    cout << "# diffractive column: full physical cross section" << endl;
    cout << "# ------------------------------------------------------------" << endl;
    cout << "# Dipole data =  " << datafile << endl;
    cout << "# sqrt(s)     =  " << param.ss << " GeV" << endl;
    cout << "# y (fixed)   =  " << param.y << endl;
    cout << "# m_charm     =  " << param.m << " GeV" << endl;
    cout << "# pt grid     =  " << n_pt << " pts in [" << pt_min << ", " << pt_max << "] GeV (log)" << endl;
    cout << "# VEGAS calls =  " << param.calls_diff << " (diffractive only; exclusive is deterministic)" << endl;
    cout << "# ============================================================" << endl;
    cout << "# K  exclusive_fixedW  diffractive_fixedW" << endl;

    vector<double> K_grid(n_pt), out_excl(n_pt), out_diff(n_pt);
    double log_lo = log(pt_min), log_hi = log(pt_max);
    for (int i = 0; i < n_pt; ++i)
        K_grid[i] = exp(log_lo + i * (log_hi - log_lo) / (n_pt - 1));

    cerr << "# Starting grid (" << n_pt << " points, "
         << omp_get_max_threads() << " OMP threads) ..." << endl;
    auto t_start = std::chrono::steady_clock::now();

    #pragma omp parallel for schedule(dynamic)
    for (int idx = 0; idx < n_pt; ++idx) {
        double K = K_grid[idx];

        parameters p = param;
        p.p  = K;
        p.p2 = K * K;
        double p2m2 = p.p2 + p.m2;
        p.denom4     = p2m2 * p2m2;
        p.p4m4       = p.p2*p.p2 + p.m2*p.m2;
        p.two_m2p2m2 = 2.0 * p.m2 * p.p2;
        p.inv_denom8 = 1.0 / (p.denom4 * p.denom4);
        p.mt   = sqrt(p2m2);

        out_excl[idx] = exclusiveCrossSection_fixedW(static_cast<void*>(&p));

        p.k_upper = p.kmax;
        out_diff[idx] = diffractiveCrossSection_fixedW(static_cast<void*>(&p));
    }

    auto t_end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();
    cerr << "# Done. Elapsed: " << (int)(elapsed/3600) << "h "
         << (int)(elapsed/60)%60 << "m " << (int)elapsed%60 << "s" << endl;

    for (int idx = 0; idx < n_pt; ++idx)
        cout << K_grid[idx] << "  " << out_excl[idx] << "  " << out_diff[idx] << endl;

    return 0;
}
