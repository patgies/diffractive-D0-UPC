#include "amplitudelib.hpp"
#include "def.hpp"
#include <gsl/gsl_errno.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

// Q_s^2 = 2/r_s^2 vs x, for N(r_s) = Ns and for the adjoint N_A = 2N - N^2.

// Adjoint Q_s^2 by bisection in ln r
static double adjoint_saturation_scale(AmplitudeLib& dipole, double x, double Ns)
{
    double lo = dipole.MinR() * 1.0001, hi = dipole.MaxR() * 0.999;
    for (int i = 0; i < 200; ++i) {
        double mid = sqrt(lo * hi);
        if (dipole.N_A(mid, x) > Ns) hi = mid; else lo = mid;
    }
    double r_s = sqrt(lo * hi);
    return 2.0 / (r_s * r_s);
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <dipole file> [Ns (default 1 - exp(-1/2))] [n_x (default 41)]" << endl;
        return 1;
    }
    gsl_set_error_handler_off();
    string datafile = argv[1];
    AmplitudeLib dipole(datafile);
    dipole.SetOutOfRangeErrors(false);

    double Ns = (argc > 2) ? atof(argv[2]) : 1.0 - exp(-0.5);
    int n_x   = (argc > 3) ? atoi(argv[3]) : 41;

    // x from XPO_MIN to XPO_FROZEN
    cout << "# Dipole data =  " << datafile << endl;
    cout << "# Q_s^2 = 2/r_s^2, N(r_s) = Ns (fundamental), N_A(r_s) = Ns (adjoint), Ns = " << Ns << endl;
    cout << "# x  Qs2_fundamental[GeV^2]  Qs2_adjoint[GeV^2]" << endl;
    cout << setprecision(6);
    double log_lo = log(XPO_MIN), log_hi = log(XPO_FROZEN);
    for (int i = 0; i < n_x; ++i) {
        // exactly XPO_FROZEN: rounding could give x > x0
        double x = (i == n_x - 1) ? XPO_FROZEN : exp(log_lo + i * (log_hi - log_lo) / (n_x - 1));
        cout << x << "  " << dipole.SaturationScale(x, Ns) << "  " << adjoint_saturation_scale(dipole, x, Ns) << endl;
    }
    return 0;
}
