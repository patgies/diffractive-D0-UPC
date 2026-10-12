#include "def.hpp"
#include "photon_flux.hpp"
#include "lhapdf_grid.hpp"
#include "bcfy_grid.hpp"
#include "kk_grid.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

// Number given by the environment variable `name`, or by `fallback` (another variable) if it is not
// set, or `value` if neither is.
static double env_number(const char* name, double value, const char* fallback = nullptr)
{
    if (getenv(name)) return atof(getenv(name));
    if (fallback && getenv(fallback)) return atof(getenv(fallback));
    return value;
}

// Settings of D0 and D0_xpom for a D0 with transverse momentum pD0 and rapidity y, read from the
// environment (see run_scripts/config.sh). d_frag_interp keeps the fragmentation function that
// param points to. Returns false, after printing the error, if a setting is not valid.
bool setup_parameters(parameters& param, unique_ptr<Interpolator>& d_frag_interp, double pD0, double y,
                      const string& frag_tag, const string& channel)
{
    load_data_and_initialize("./input/WS_photon_flux/Gamma_AA.dat");

    // TARGET=pA: proton target. It changes the collision energy and the photon flux.
    param.target = getenv("TARGET") ? getenv("TARGET") : "AA";
    param.ss     = (param.target == "pA") ? 8160.0 : 5360.0;

    // VEGAS calls per point. Exclusive: 5D integral below EXCL_PT_THRESHOLD, 3D integral above.
    param.calls_excl            = (size_t)env_number("CALLS_EXCL", 1e5, "CALLS");
    param.calls_excl_factorized = (size_t)env_number("CALLS_EXCL_FACTORIZED", 2e3);
    param.excl_pt_threshold     = env_number("EXCL_PT_THRESHOLD", 3.0);
    param.calls_diff            = (size_t)env_number("CALLS_DIFF", 1e5, "CALLS");
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
            return false;
        }
        init_pA_flux(env_number("SIGMA_NN", 99.0));
        if (param.flux_model == "WS") init_ws_form_factor();
        const double R_PL_fm = 7.1;   // point-like radius, as in b_min = 2R_PL = 14.2 fm for Pb+Pb (arXiv:2606.05469)
        param.bmin = (param.flux_model == "PL") ? (1.1 * R_PL_fm / 0.197327) : 1e-3;
    } else {
        param.flux_model = getenv("FLUX_MODEL") ? getenv("FLUX_MODEL") : "EFF";
        if      (param.flux_model == "WS")  init_ws_form_factor();
        else if (param.flux_model == "EFF") init_effective_flux(param.channel, &param);
        else if (param.flux_model == "STARLIGHT") init_starlight_flux(param.channel);
    }
    // QPMAX: largest photon energy, only for tests (default 800 GeV).
    param.qpmax   = env_number("QPMAX", 800.0);

    // Upper limit of the b integral, 60/(z_gamma_min*mn), where z_gamma_min is the smallest
    // photon energy fraction possible at this (pD0, y).
    double mt_min = sqrt(pD0*pD0 + param.m2);          // mt at z_h = 1
    double z_gamma_min = mt_min * exp(y) / param.ss;
    param.bmax    = 60.0 / (z_gamma_min * param.mn);

    // Fragmentation (c -> D0)
    param.z_h_min = 0.05;
    param.z_h_max = 1.0;

    if (frag_tag == "BCFY")              param.frag_type = FragmentationType::BCFY;
    else if (frag_tag == "KniehlKramer") param.frag_type = FragmentationType::KniehlKramer;
    else if (frag_tag == "HymnD")        param.frag_type = FragmentationType::HymnD;
    else {
        cerr << "Error: unknown frag_type '" << frag_tag << "'. Expected BCFY, KniehlKramer or HymnD." << endl;
        return false;
    }

    string hymnD_file = getenv("HYMND_FILE") ? getenv("HYMND_FILE") : "input/HymnD/prompt-D0-1-109_0000.dat";
    const int hymnD_charm_flavor = 4;  // PDG id for the charm quark
    // Fragmentation scale Q = SCALE_FACTOR * mt0, mt0 = sqrt(pD0^2 + m^2).
    double mt0 = sqrt(pD0*pD0 + param.m2);
    double frag_scale = env_number("SCALE_FACTOR", 1.0) * mt0;
    // the fragmentation scale cannot be below the charm mass
    if (frag_scale < param.m) frag_scale = param.m;
    if (param.frag_type == FragmentationType::HymnD) {
        d_frag_interp = MakeLHAPDFGridInterpolator(hymnD_file, hymnD_charm_flavor, frag_scale);
    } else if (param.frag_type == FragmentationType::BCFY) {
        d_frag_interp = MakeBCFYInterpolator(frag_scale);
    } else {
        d_frag_interp = MakeKniehlKramerInterpolator(frag_scale);
    }
    param.D_frag_interp = d_frag_interp.get();

    param.pD0 = pD0;
    param.y   = y;
    return true;
}
