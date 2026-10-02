#include "amplitudelib.hpp"
#include "def.hpp"
#include "photon_flux.hpp"
#include "fragmentation.hpp"
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


//  Run: ./D0_xpom <dipole_file> <pD0> <y> <x_po> [frag_type] [channel]
//   x_po       pomeron momentum fraction, should be between 0 and 0.1
//   frag_type  BCFY (default) | KniehlKramer | HymnD
//   channel    An0n (default) | Xn0n | 0n0n | PL(AnAn)
//
// This is the x_po-fixed version of D0, giving both channels differential
// in x_P instead of integrated over it. If you run this for a bunch of
// x_po values and integrate the results over x_po (e.g. with Simpson's
// rule), you should get back the same numbers D0 gives you for both
// columns. For diffractive, x_po is a genuinely independent variable (still
// held fixed while q+ is separately integrated); for exclusive, x_P is not
// independent -- it's a function of q+ (and y), so "fixed x_P" instead
// picks out the one q+ that gives it and applies the resulting Jacobian
// (see integrand_exclusive_xpom in integrand.cpp).


int main(int argc, char* argv[])
{
    if (argc < 5 || argc > 7) {
        cerr << "Usage: " << argv[0] << " <dipole_file> <pD0> <y> <x_po> [frag_type] [channel]" << endl;
        return 1;
    }
    string datafile  = argv[1];
    double pD0       = atof(argv[2]);
    double y         = atof(argv[3]);
    double x_po      = atof(argv[4]);
    string frag_tag;
    if (argc > 5) frag_tag = argv[5];
    else          frag_tag = "BCFY";

    string channel;
    if (argc > 6) channel = argv[6];
    else          channel = "An0n";

    AmplitudeLib inst(datafile);
    inst.SetOutOfRangeErrors(false);
    inst.SetInterpolationMethod(LINEAR_LINEAR);  

    gsl_set_error_handler_off();
    load_data_and_initialize("./inputs/Gamma_AA.dat");

    parameters param;
    param.dipole   = &inst;
    param.datafile = datafile;
    // TARGET=pA (proton target, Sec. 6 of arXiv:2606.05469) -- see main.cpp
    // for the full explanation.
    param.target = getenv("TARGET") ? getenv("TARGET") : "AA";
    param.ss     = (param.target == "pA") ? 8160.0 : 5360.0;

    const char* calls_default = getenv("CALLS");
    const char* calls_excl_fact = getenv("CALLS_EXCL_FACTORIZED");
    const char* calls_diff    = getenv("CALLS_DIFF");
    param.calls_excl = (size_t)atof(calls_default ? calls_default : "1e5");
    param.calls_diff = (size_t)atof(calls_diff ? calls_diff : (calls_default ? calls_default : "1e5"));
    // exclusiveCrossSection_xpom uses the factorized (analytic r1,r2)
    // integrand exclusively, same as exclusiveCrossSection above threshold
    // -- see main.cpp for the same default.
    param.calls_excl_factorized = (size_t)atof(calls_excl_fact ? calls_excl_fact : "1e3");
    param.m        = 1.5;
    param.m2       = param.m * param.m;

    param.alpha   = 1.0/137.0;
    param.Z       = 82.0;
    param.mn      = 0.938;
    param.S       = pow(17.4, 2) / pow(0.197327, 2);
    param.channel = channel;
    param.bmin    = 14.2 / 0.197327;
    // FLUX_MODEL: "EFF" (default) is the effective-flux geometric
    // convolution; "PL" is the old single-b treatment; "WS" is identical to
    // "PL" here -- see main.cpp for the full explanation. TARGET=pA: see
    // main.cpp for the pA flux setup (default flux_model there is "WS").
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
    param.qpmax   = 800.0;

    // Adaptive bmax, replacing the old fixed 650 -- see main.cpp for the
    // derivation (same (pD0, y) kinematics here, so the same worst-case
    // smallest-z argument applies).
    double mt_min = sqrt(pD0*pD0 + param.m2);
    double z_min  = mt_min * exp(y) / param.ss;
    param.bmax    = 60.0 / (z_min * param.mn);

    // Fragmentation (c -> D0)
    param.r      = 0.1;
    param.N_kk   = 0.694;
    param.eps_kk = 0.101;
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
    else                       hymnD_file = "inputs/prompt-D0-1-109/prompt-D0-1-109_0000.dat";
    const int hymnD_charm_flavor = 4;  // PDG id for the charm quark
    // Scale-variation knob: fragmentation scale Q = scale_factor * mt0,
    // mt0 = sqrt(pD0^2 + m^2) the D0 transverse mass (see main.cpp).
    double mt0 = sqrt(pD0*pD0 + param.m2);
    double scale_factor = getenv("SCALE_FACTOR") ? atof(getenv("SCALE_FACTOR")) : 1.0;
    double frag_scale = scale_factor * mt0;
    // Below the charm mass, the fragmentation function is undefined (the
    // HymnD grid's charm production threshold): use the charm mass itself
    // as a floor rather than letting a small scale_factor push Q below it.
    if (frag_scale < param.m) frag_scale = param.m;
    // All three frag_types are DGLAP-evolved z-interpolators at frag_scale
    // now (see main.cpp for the same convention).
    unique_ptr<Interpolator> d_frag_interp;
    if (param.frag_type == FragmentationType::HymnD) {
        d_frag_interp = MakeHymnDZInterpolator(hymnD_file, hymnD_charm_flavor, frag_scale);
    } else if (param.frag_type == FragmentationType::BCFY) {
        d_frag_interp = MakeBCFYInterpolator(frag_scale);
    } else {
        d_frag_interp = MakeKniehlKramerInterpolator(frag_scale);
    }
    param.D_frag_interp = d_frag_interp.get();

    param.pD0       = pD0;
    param.y         = y;
    param.fixed_xpo = x_po;

    cout << "# fragmentation : " << frag_tag << endl;
    cout << "# target        : " << param.target << " (sqrt(s_NN)=" << param.ss << " GeV, flux_model=" << param.flux_model << ")" << endl;
    cout << "# channel       : " << param.channel << (param.target == "pA" ? " (unused: no EMD for pA)" : "") << endl;
    cout << "# dipole file   : " << datafile << endl;
    cout << "# pD0  x_po  exclusive_dxpo  diffractive_dxpo" << endl;

    double result_excl = exclusiveCrossSection_xpom(static_cast<void*>(&param));
    double result_diff = diffractiveCrossSection_xpom(static_cast<void*>(&param));

    cout << pD0 << "  " << x_po << "  " << result_excl << "  " << result_diff << endl;

    return 0;
}
