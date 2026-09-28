#include "amplitudelib.hpp"
#include "def.hpp"
#include "photon_flux.hpp"
#include <string>
#include <iostream>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>

using namespace std;

// Fixed-W (no photon-flux), x_po-differential diffractive spectrum, as a
// function of x_po, at a fixed rapidity y and observed quark transverse
// momentum K. Same "q+ = 3p+ fixed by hand" trick as main_fixedW.cpp, but
// x_po is fixed per point (like main_xpom.cpp) instead of being integrated
// over -- the fixed-W analogue of D0_xpom.
//
// Run: ./D0_fixedW_xpom <y> <K> <x_po> [calls]
//
// Exclusive doesn't depend on x_po (no rapidity-gap variable in that
// process), so it's printed unchanged from exclusiveCrossSection_fixedW
// for reference, same "missing prefactor" convention as main_fixedW.cpp.

int main(int argc, char* argv[])
{
    if (argc < 4 || argc > 5) {
        cerr << "Usage: " << argv[0] << " <y> <K> <x_po> [calls]" << endl;
        return 1;
    }

    string datafile = "./data/proton/mve.dat";

    AmplitudeLib inst(datafile);
    inst.SetOutOfRangeErrors(false);

    gsl_set_error_handler_off();
    gsl_rng_env_setup();

    load_data_and_initialize("./inputs/Gamma_AA.dat");   // unused (no flux here)

    parameters param;
    param.dipole   = &inst;
    param.datafile = datafile;
    param.ss       = 5360.0;
    param.calls_diff = (argc > 4) ? (size_t)atof(argv[4]) : (size_t)1e5;
    param.calls_excl = param.calls_diff;
    param.m        = 1.5;
    param.m2       = param.m * param.m;
    param.y        = atof(argv[1]);
    param.kmax     = 400.0;

    double K = atof(argv[2]);
    param.fixed_xpo = atof(argv[3]);

    param.p  = K;
    param.p2 = K * K;
    double p2m2 = param.p2 + param.m2;
    param.denom4     = p2m2 * p2m2;
    param.p4m4       = param.p2*param.p2 + param.m2*param.m2;
    param.two_m2p2m2 = 2.0 * param.m2 * param.p2;
    param.inv_denom8 = 1.0 / (param.denom4 * param.denom4);
    param.mt   = sqrt(p2m2);
    param.k_upper = param.kmax;

    double result_excl = exclusiveCrossSection_fixedW(static_cast<void*>(&param));
    double result_diff = diffractiveCrossSection_fixedW_xpom(static_cast<void*>(&param));

    cout << "# Fixed-W (no photon flux), x_po-differential diffractive dsigma/(d2K dx_po)" << endl;
    cout << "# y (fixed)  =  " << param.y << endl;
    cout << "# K (fixed)  =  " << K << endl;
    cout << "# exclusive column:   missing prefactor alpha*Nc*ef^2*Sperp/(2pi^2), independent of x_po" << endl;
    cout << "# diffractive column: full physical cross section at this x_po" << endl;
    cout << "# x_po  exclusive_fixedW  diffractive_fixedW_dxpo" << endl;
    cout << param.fixed_xpo << "  " << result_excl << "  " << result_diff << endl;

    return 0;
}
