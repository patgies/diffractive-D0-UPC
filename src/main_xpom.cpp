#include "amplitudelib.hpp"
#include "def.hpp"
#include <string>
#include <iostream>
#include <cstdlib>
#include <memory>
#include <gsl/gsl_errno.h>

using namespace std;

// As D0, at fixed x_po: exclusive dsigma/dx_po, diffractive dsigma/dln(x_po).
int main(int argc, char* argv[])
{
    if (argc < 5 || argc > 7) {
        cerr << "Usage: " << argv[0] << " <dipole_file> <pD0> <y> <x_po> [frag_type] [channel]" << endl;
        return 1;
    }
    string datafile = argv[1];
    double pD0      = atof(argv[2]);
    double y        = atof(argv[3]);
    double x_po     = atof(argv[4]);
    string frag_tag = (argc > 5) ? argv[5] : "BCFY";
    string channel  = (argc > 6) ? argv[6] : "An0n";

    // GSL errors off (see main.cpp)
    gsl_set_error_handler_off();
    AmplitudeLib inst(datafile);
    inst.SetOutOfRangeErrors(false);
    inst.SetInterpolationMethod(LINEAR_LINEAR);

    if (getenv("DIPOLE_X0")) {
        inst.SetX0(atof(getenv("DIPOLE_X0")));
    }

    parameters param;
    param.dipole = &inst;
    unique_ptr<Interpolator> d_frag_interp;
    if (!setup_parameters(param, d_frag_interp, pD0, y, frag_tag, channel)) return 1;
    param.fixed_xpo = x_po;

    cout << "# fragmentation : " << frag_tag << endl;
    cout << "# target        : " << param.target << " (sqrt(s_NN)=" << param.ss << " GeV, flux_model=" << param.flux_model << ")" << endl;
    cout << "# channel       : " << param.channel << (param.target == "pA" ? " (unused: no EMD for pA)" : "") << endl;
    cout << "# dipole file   : " << datafile << endl;
    cout << "# DTMD          : " << param.dtmd << " (gluon, quark or both)" << endl;
    cout << "# pD0  x_po  exclusive_dxpo  diffractive_dlnxpo" << endl;

    double result_excl = exclusiveCrossSection_xpom(static_cast<void*>(&param));
    double result_diff = diffractiveCrossSection_xpom(static_cast<void*>(&param));

    cout << pD0 << "  " << x_po << "  " << result_excl << "  " << result_diff << endl;

    return 0;
}
