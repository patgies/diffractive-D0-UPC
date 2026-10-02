#include "amplitudelib.hpp"
#include "def.hpp"
#include "photon_flux.hpp"
#include "hymnd_grid.hpp"
#include "bcfy_grid.hpp"
#include "kk_grid.hpp"
#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <gsl/gsl_errno.h>

using namespace std;


// Run: ./D0 <dipole_file> <pD0> <y> [frag_type] [channel]
// frag_type: BCFY | KniehlKramer | HymnD;  channel: An0n | Xn0n | 0n0n | PL(AnAn)


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

    // DIPOLE_X0 sets x0 for dipole files that do not have it in the header.
    if (getenv("DIPOLE_X0")) {
        inst.SetX0(atof(getenv("DIPOLE_X0")));
    }

    gsl_set_error_handler_off();
    load_data_and_initialize("./inputs/WS_photon_flux/Gamma_AA.dat");

    parameters param;
    param.dipole   = &inst;
    // TARGET=pA: proton target. It changes the collision energy and the photon flux.
    param.target = getenv("TARGET") ? getenv("TARGET") : "AA";
    param.ss     = (param.target == "pA") ? 8160.0 : 5360.0;


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
    // FLUX_MODEL: EFF (default, effective flux), PL or WS (flux at one b),
    // STARLIGHT (from a table). For TARGET=pA: WS (default) or PL.
    if (param.target == "pA") {
        param.flux_model = getenv("FLUX_MODEL") ? getenv("FLUX_MODEL") : "WS";
        if (param.flux_model == "EFF" || param.flux_model == "STARLIGHT") {
            cerr << "Error: FLUX_MODEL=" << param.flux_model << " is not meaningful with TARGET=pA "
                    "(no target-nucleus extent to convolve over). Use PL or WS." << endl;
            return 1;
        }
        double sigma_NN_mb = getenv("SIGMA_NN") ? atof(getenv("SIGMA_NN")) : 99.0;
        init_pA_flux(sigma_NN_mb);
        if (param.flux_model == "WS") init_ws_form_factor();
        const double RA_fm = 6.49;
        param.bmin = (param.flux_model == "PL") ? (1.1 * RA_fm / 0.197327) : 1e-3;
    } else {
        param.flux_model = getenv("FLUX_MODEL") ? getenv("FLUX_MODEL") : "EFF";
        if      (param.flux_model == "WS")  init_ws_form_factor();
        else if (param.flux_model == "EFF") init_effective_flux(param.channel, &param);
        else if (param.flux_model == "STARLIGHT") init_starlight_flux(param.channel);
    }
    // QPMAX: largest photon energy, only for tests (default 800 GeV).
    param.qpmax   = getenv("QPMAX") ? atof(getenv("QPMAX")) : 800.0;

    // Upper limit of the b integral, 60/(z_min*mn), where z_min is the smallest
    // photon energy fraction possible at this (pD0, y).
    double mt_min = sqrt(pD0*pD0 + param.m2);       // mt at zh=zmax=1
    double z_min  = mt_min * exp(y) / param.ss;     // smallest photon energy fraction possible
    param.bmax    = 60.0 / (z_min * param.mn);

    // Fragmentation (c -> D0)
    param.zmin   = 0.05;
    param.zmax   = 1.0;

    if (frag_tag == "BCFY")              param.frag_type = FragmentationType::BCFY;
    else if (frag_tag == "KniehlKramer") param.frag_type = FragmentationType::KniehlKramer;
    else if (frag_tag == "HymnD")       param.frag_type = FragmentationType::HymnD;
    else {
        cerr << "Error: unknown frag_type '" << frag_tag << "'. Expected BCFY, KniehlKramer or HymnD." << endl;
        return 1;
    }

    string hymnD_file;
    if (getenv("HYMND_FILE")) hymnD_file = getenv("HYMND_FILE");
    else                       hymnD_file = "inputs/HymnD/prompt-D0-1-109_0000.dat";
    const int hymnD_charm_flavor = 4;  // PDG id for the charm quark
    // Fragmentation scale Q = SCALE_FACTOR * mt0, mt0 = sqrt(pD0^2 + m^2).
    double mt0 = sqrt(pD0*pD0 + param.m2);
    double scale_factor = getenv("SCALE_FACTOR") ? atof(getenv("SCALE_FACTOR")) : 1.0;
    double frag_scale = scale_factor * mt0;
    // the fragmentation scale cannot be below the charm mass
    if (frag_scale < param.m) frag_scale = param.m;
    // D(z) at the fragmentation scale for the chosen frag_type.
    unique_ptr<Interpolator> d_frag_interp;
    if (param.frag_type == FragmentationType::HymnD) {
        d_frag_interp = MakeHymnDZInterpolator(hymnD_file, hymnD_charm_flavor, frag_scale);
    } else if (param.frag_type == FragmentationType::BCFY) {
        d_frag_interp = MakeBCFYInterpolator(frag_scale);
    } else {
        d_frag_interp = MakeKniehlKramerInterpolator(frag_scale);
    }
    param.D_frag_interp = d_frag_interp.get();

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
    cout << "# target        : " << param.target << " (sqrt(s_NN)=" << param.ss << " GeV, flux_model=" << param.flux_model << ")" << endl;
    cout << "# channel       : " << param.channel << (param.target == "pA" ? " (unused: no EMD for pA)" : "") << endl;
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
