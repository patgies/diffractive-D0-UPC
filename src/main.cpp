#include "amplitudelib.hpp"
#include "def.hpp"
#include "gamma_aa.hpp"
#include "fragmentation.hpp"
#include "lhapdf_grid.hpp"
#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <gsl/gsl_errno.h>

using namespace std;


//  Run: ./D0 <dipole_file> <pD0> <y> [frag_type] [channel]
//   frag_type  BCFY (default) | KniehlKramer | LHAPDF
//   channel    An0n (default) | Xn0n | PL(AnAn)


int main(int argc, char* argv[])
{
    if (argc < 4 || argc > 6) {
        cerr << "Usage: " << argv[0] << " <dipole_file> <pD0> <y> [frag_type] [channel]" << endl;
        return 1;
    }
    string datafile  = argv[1];
    double pD0       = atof(argv[2]);
    double y         = atof(argv[3]);
    string frag_tag;
    if (argc > 4) frag_tag = argv[4];
    else          frag_tag = "BCFY";

    string channel;
    if (argc > 5) channel = argv[5];
    else          channel = "An0n";

    AmplitudeLib inst(datafile);
    inst.SetOutOfRangeErrors(false);
    inst.SetInterpolationMethod(LINEAR_LINEAR);   

    gsl_set_error_handler_off();
    load_data_and_initialize("./data/Gamma_AA.dat");

    parameters param;
    param.dipole   = &inst;
    param.datafile = datafile;
    param.ss       = 5360.0;

 
    const char* calls_default = getenv("CALLS");
    const char* calls_excl    = getenv("CALLS_EXCL");
    const char* calls_excl_fact = getenv("CALLS_EXCL_FACTORIZED");
    const char* calls_diff    = getenv("CALLS_DIFF");
    const char* excl_pt_threshold = getenv("EXCL_PT_THRESHOLD");
    if (calls_excl)         param.calls_excl = (size_t)atof(calls_excl);
    else if (calls_default) param.calls_excl = (size_t)atof(calls_default);
    else                     param.calls_excl = (size_t)1e5;

    if (calls_excl_fact) param.calls_excl_factorized = (size_t)atof(calls_excl_fact);
    else                  param.calls_excl_factorized = (size_t)1e3;

    if (excl_pt_threshold) param.excl_pt_threshold = atof(excl_pt_threshold);
    else                    param.excl_pt_threshold = 4.0;

    if (calls_diff)         param.calls_diff = (size_t)atof(calls_diff);
    else if (calls_default) param.calls_diff = (size_t)atof(calls_default);
    else                     param.calls_diff = (size_t)1e5;
    param.m        = 1.5;
    param.m2       = param.m * param.m;


    param.alpha   = 1.0/137.0;
    param.Z       = 82.0;
    param.mn      = 0.938;
    param.S       = pow(17.4, 2) / pow(0.197327, 2);
    param.channel = channel;
    param.bmin    = 14.2 / 0.197327;
    param.bmax    = 650.0;
    param.qpmax   = 800.0;

    // Fragmentation (c -> D0)
    param.r      = 0.1;
    param.N_kk   = 0.694;
    param.eps_kk = 0.101;
    param.zmin   = 0.05;
    param.zmax   = 1.0;

    if (frag_tag == "BCFY")              param.frag_type = FragmentationType::BCFY;
    else if (frag_tag == "KniehlKramer") param.frag_type = FragmentationType::KniehlKramer;
    else if (frag_tag == "LHAPDF")       param.frag_type = FragmentationType::LHAPDF;
    else {
        cerr << "Error: unknown frag_type '" << frag_tag << "'. Expected BCFY, KniehlKramer or LHAPDF." << endl;
        return 1;
    }

    string lhapdf_file;
    if (getenv("LHAPDF_FILE")) lhapdf_file = getenv("LHAPDF_FILE");
    else                       lhapdf_file = "data/prompt-D0-1-109/prompt-D0-1-109_0000.dat";
    const int lhapdf_charm_flavor = 4;  // PDG id for the charm quark
    unique_ptr<Interpolator> d_frag_interp;  
    if (param.frag_type == FragmentationType::LHAPDF) {
        d_frag_interp = MakeLHAPDFZInterpolator(lhapdf_file, lhapdf_charm_flavor, param.m);
        param.D_frag_interp = d_frag_interp.get();
    }

    param.pD0 = pD0;
    param.y   = y;



    const char* process_env = getenv("PROCESS");
    string process_mode;
    if (process_env) process_mode = process_env;
    else              process_mode = "both";
    if (process_mode != "both" && process_mode != "exclusive" && process_mode != "diffractive") {
        cerr << "Error: unknown PROCESS '" << process_mode << "'. Expected both, exclusive or diffractive." << endl;
        return 1;
    }

    cout << "# fragmentation : " << frag_tag << endl;
    cout << "# channel       : " << param.channel << endl;
    cout << "# dipole file   : " << datafile << endl;
    cout << "# process       : " << process_mode << endl;
    cout << "# pD0  exclusive  diffractive" << endl;

    double result_excl = 0.0, result_diff = 0.0;
    if (process_mode == "both" || process_mode == "exclusive")
        result_excl = exclusiveCrossSection(static_cast<void*>(&param));
    if (process_mode == "both" || process_mode == "diffractive")
        result_diff = diffractiveCrossSection(static_cast<void*>(&param));

    cout << pD0 << "  " << result_excl << "  " << result_diff << endl;

    return 0;
}
