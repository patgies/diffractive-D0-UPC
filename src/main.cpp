#include "amplitudelib.hpp"
#include "def.hpp"
#include <string>
#include <iostream>
#include <cstdlib>
#include <memory>
#include <gsl/gsl_errno.h>

using namespace std;

// D0 cross section at one (pD0, y), without prefactors. Settings: run_scripts/config.sh, setup.cpp.
int main(int argc, char* argv[])
{
    if (argc < 4 || argc > 6) {
        cerr << "Usage: " << argv[0] << " <dipole_file> <pD0> <y> [frag_type] [channel]" << endl;
        return 1;
    }
    string datafile = argv[1];
    double pD0      = atof(argv[2]);
    double y        = atof(argv[3]);
    string frag_tag = (argc > 4) ? argv[4] : "BCFY";
    string channel  = (argc > 5) ? argv[5] : "An0n";

    // GSL errors off: the dilute-edge files have no saturation scale
    gsl_set_error_handler_off();
    AmplitudeLib inst(datafile);
    inst.SetOutOfRangeErrors(false);
    inst.SetInterpolationMethod(LINEAR_LINEAR);

    // DIPOLE_X0: x0 if not in the file header
    if (getenv("DIPOLE_X0")) {
        inst.SetX0(atof(getenv("DIPOLE_X0")));
    }

    parameters param;
    param.dipole = &inst;
    unique_ptr<Interpolator> d_frag_interp;
    if (!setup_parameters(param, d_frag_interp, pD0, y, frag_tag, channel)) return 1;

    string process_mode = getenv("PROCESS") ? getenv("PROCESS") : "both";
    if (process_mode != "both" && process_mode != "exclusive" && process_mode != "diffractive") {
        cerr << "Error: unknown PROCESS '" << process_mode << "'. Expected both, exclusive or diffractive." << endl;
        return 1;
    }

    cout << "# fragmentation : " << frag_tag << endl;
    cout << "# target        : " << param.target << " (sqrt(s_NN)=" << param.ss << " GeV, flux_model=" << param.flux_model << ")" << endl;
    cout << "# channel       : " << param.channel << (param.target == "pA" ? " (unused: no EMD for pA)" : "") << endl;
    cout << "# dipole file   : " << datafile << endl;
    cout << "# process       : " << process_mode << endl;
    cout << "# DTMD          : " << param.dtmd << " (gluon, quark or both)" << endl;
    cout << "# pD0  exclusive  diffractive" << endl;

    double result_excl = 0.0, result_diff = 0.0;
    if (process_mode == "both" || process_mode == "exclusive")
        result_excl = exclusiveCrossSection(static_cast<void*>(&param));
    if (process_mode == "both" || process_mode == "diffractive")
        result_diff = diffractiveCrossSection(static_cast<void*>(&param));

    cout << pD0 << "  " << result_excl << "  " << result_diff << endl;

    return 0;
}
