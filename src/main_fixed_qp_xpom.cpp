#include "amplitudelib.hpp"
#include "def.hpp"
#include <cstdlib>
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

// Diffractive dsigma/dln(x_po) at fixed q+ (qp = 2*pp, no photon flux), at fixed rapidity y and
// quark transverse momentum p, on a log grid in x_po. The exclusive part has no x_po
// distribution here: at fixed q+ it sits at x_po = M^2/W^2, where the diffractive one starts.

int main(int argc, char* argv[])
{
    // arguments: y, p, number of x_po points, VEGAS calls, dipole file (default: proton)
    string datafile = (argc > 5) ? argv[5] : "./data/proton/mve.dat";

    // GSL errors off before loading the dipole (see main.cpp)
    gsl_set_error_handler_off();
    AmplitudeLib inst(datafile);
    inst.SetOutOfRangeErrors(false);

    if (getenv("DIPOLE_X0")) {
        inst.SetX0(atof(getenv("DIPOLE_X0")));
    }

    gsl_rng_env_setup();

    parameters param;
    param.dipole   = &inst;
    param.ss       = 5360.0;
    param.calls_diff = (argc > 4) ? (size_t)atof(argv[4]) : (size_t)1e5;
    param.m        = 1.5;
    param.m2       = param.m * param.m;
    param.y        = (argc > 1) ? atof(argv[1]) : 0.0;
    param.kmax     = 400.0;

    param.p    = (argc > 2) ? atof(argv[2]) : 2.0;
    param.p2   = param.p * param.p;
    double p2m2 = param.p2 + param.m2;
    param.mt   = sqrt(p2m2);
    param.k_upper = param.kmax;

    int n_xpo = (argc > 3) ? atoi(argv[3]) : 40;
    double xpo_min = 1e-5, xpo_max = 0.1;   // x_po > 0.1 is zero, as in charm_fixed_qp

    // smallest x_po of the diffractive part (x < 1): M^2/W^2 with z = 1/2
    double xpo_excl = 4.0 * p2m2 / (M_SQRT2 * 2.0 * (param.mt / M_SQRT2) * exp(param.y) * param.ss);

    std::time_t now = std::time(nullptr);
    char date_buf[32];
    std::strftime(date_buf, sizeof(date_buf), "%Y-%m-%d %H:%M:%S", std::localtime(&now));

    cout << "# ============================================================" << endl;
    cout << "# Date        =  " << date_buf << endl;
    cout << "# Fixed-q+ (no photon flux) diffractive dsigma/(d2p dln x_po)" << endl;
    cout << "# diffractive   : gluon dipole N_g, per d ln(x_po)" << endl;
    cout << "# W^2 = sqrt(2)*ss*2*pp fixed via qp=2*pp (z=1/2); y and p fixed" << endl;
    cout << "# diffractive column: missing prefactor alpha_s*alpha*ef^2*(Nc^2-1)*Sperp/(8pi^4)" << endl;
    cout << "# ------------------------------------------------------------" << endl;
    cout << "# Dipole data =  " << datafile << endl;
    cout << "# sqrt(s)     =  " << param.ss << " GeV" << endl;
    cout << "# y (fixed)   =  " << param.y << endl;
    cout << "# p (fixed)   =  " << param.p << endl;
    cout << "# m_charm     =  " << param.m << " GeV" << endl;
    cout << "# x_po (excl) =  " << xpo_excl << " (exclusive x_po, start of the diffractive one)" << endl;
    cout << "# x_po grid   =  " << n_xpo << " pts in [" << xpo_min << ", " << xpo_max << "] (log)" << endl;
    cout << "# VEGAS calls =  " << param.calls_diff << endl;
    cout << "# ============================================================" << endl;
    cout << "# x_po  diffractive_fixed_qp_dlnxpo" << endl;

    vector<double> xpo_grid(n_xpo), out_diff(n_xpo);
    double log_lo = log(xpo_min), log_hi = log(xpo_max);
    for (int i = 0; i < n_xpo; ++i)
        xpo_grid[i] = exp(log_lo + i * (log_hi - log_lo) / (n_xpo - 1));

    cerr << "# Starting grid (" << n_xpo << " points, "
         << omp_get_max_threads() << " OMP threads) ..." << endl;
    auto t_start = std::chrono::steady_clock::now();

    #pragma omp parallel for schedule(dynamic)
    for (int idx = 0; idx < n_xpo; ++idx) {
        parameters p = param;
        p.fixed_xpo = xpo_grid[idx];
        out_diff[idx] = (p.fixed_xpo <= xpo_excl) ? 0.0
                        : diffractiveCrossSection_fixed_qp_xpom(static_cast<void*>(&p));
    }

    auto t_end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();
    cerr << "# Done. Elapsed: " << (int)(elapsed/3600) << "h "
         << (int)(elapsed/60)%60 << "m " << (int)elapsed%60 << "s" << endl;

    for (int idx = 0; idx < n_xpo; ++idx)
        cout << xpo_grid[idx] << "  " << out_diff[idx] << endl;

    return 0;
}
