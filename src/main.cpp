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


//  Run: ./D0 <dipole_file> <pD0> <y> [frag_type] [channel]
//   frag_type  BCFY (default) | KniehlKramer | HymnD
//   channel    An0n (default) | Xn0n | 0n0n | PL(AnAn)


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

    // Some BK solver output (e.g. the rcbk-produced BK-initial-condition
    // posterior-sample dipole files under bk/) leaves x0 out of the header,
    // which DataFile reads as an "invalid x0" and resets to 0 -- fix that up
    // here rather than editing the data files. Only takes effect if
    // explicitly requested, so it's a no-op for every other dataset's
    // already-correct header value.
    if (getenv("DIPOLE_X0")) {
        inst.SetX0(atof(getenv("DIPOLE_X0")));
    }

    gsl_set_error_handler_off();
    load_data_and_initialize("./inputs/Gamma_AA.dat");

    parameters param;
    param.dipole   = &inst;
    param.datafile = datafile;
    // TARGET=pA (proton target, Sec. 6 of arXiv:2606.05469) switches both
    // the collision energy and the photon-flux geometry; see the flux_model
    // block below for the rest of the pA setup.
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
    // FLUX_MODEL: "EFF" (default) is the target-nucleus geometric
    // convolution (arXiv:2404.09731 Eq. 4, effective_photon_flux() in
    // photon_flux.cpp) -- validated in scan_flux_eff.cpp and
    // photon_flux_discrepancy.pdf to actually reproduce Paakkinen's
    // tabulated Starlight flux, unlike the old single-b treatment. "PL" is
    // that old single-b treatment (Eq. 1/15 of our own paper and of Guzey
    // et al. 2606.05469), kept for comparison/reproducing old results;
    // "WS" swaps in the Woods-Saxon bare flux (flux_density_WS()) within
    // that same single-b treatment, but is numerically IDENTICAL to "PL"
    // there, since par.bmin always exceeds R_A (see flux_density_WS()'s
    // own comment) -- it's "EFF" that actually fixes the bmin-independent
    // part of the discrepancy, not "WS".
    // TARGET=pA: EFF doesn't apply (Sec. 6.1 of arXiv:2606.05469 uses the
    // single-b treatment, since a proton target has no extent of its own to
    // convolve over), so the default there is "WS" instead, with
    // Gamma_pA(b)=exp(-sigma_NN*T_A(b)) (init_pA_flux()/GammaPA()) standing
    // in for Gamma_AA -- see photon_flux()'s target=="pA" branch. "PL" is
    // the sharp point-like comparison cutoff at bmin=1.1*R_A (their Fig. 11).
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
    // QPMAX override: diagnostic only, to check what fraction of the cross
    // section comes from photon energies above a given y=2qp/sqrt(2)/ss
    // threshold (compare a run with QPMAX=800 against one capped at the
    // qp corresponding to that y). Defaults to the physical 800 GeV.
    param.qpmax   = getenv("QPMAX") ? atof(getenv("QPMAX")) : 800.0;

    // Adaptive upper limit for the b integral (replaces a fixed bmax=650):
    // large enough that the photon flux (photon_flux()/flux_density() in
    // photon_flux.cpp, which already dies off past eta=z*mn*b=50) has fully
    // decayed before the box edge, for every point VEGAS can sample here --
    // same idea as the bmax=60/(z*mn) used to validate the flux against
    // arXiv:2404.09731 (src/scan_flux.cpp), but z there is the photon
    // energy fraction, NOT this D0's rapidity y, so it can't be plugged in
    // directly: q+ (and thus z) is itself a Monte Carlo variable sampled
    // inside the integrand (integrand.cpp), not fixed externally like in
    // scan_flux.cpp. Instead we use the smallest z reachable at this
    // (pD0, y): q+ >= p+ = (mt/sqrt2)*exp(y), and mt=sqrt(pc^2+m^2) with
    // pc=pD0/zh is smallest at zh=zmax=1 (pc=pD0), so that fixes the
    // worst-case (smallest) z this run will ever sample.
    double mt_min = sqrt(pD0*pD0 + param.m2);       // mt at zh=zmax=1
    double z_min  = mt_min * exp(y) / param.ss;     // smallest photon energy fraction reachable
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
    // mt0 = sqrt(pD0^2 + m^2) the D0 transverse mass (matching
    // inclusive-D0-UPC's convention, not just the fixed charm mass).
    // Defaults to scale_factor=1.0 (today's behavior). Set to 0.5/2.0 for
    // the conventional up/down envelope around that central scale -- a
    // separate uncertainty source from the HymnD replica band, so it only
    // needs the central member (0000), not all 101 replicas.
    double mt0 = sqrt(pD0*pD0 + param.m2);
    double scale_factor = getenv("SCALE_FACTOR") ? atof(getenv("SCALE_FACTOR")) : 1.0;
    double frag_scale = scale_factor * mt0;
    // Below the charm mass, the fragmentation function is undefined (the
    // HymnD grid's charm production threshold): use the charm mass itself
    // as a floor rather than letting a small scale_factor push Q below it.
    if (frag_scale < param.m) frag_scale = param.m;
    // All three frag_types are DGLAP-evolved z-interpolators at frag_scale
    // now: BCFY/KniehlKramer from their own eko-evolved grids (see
    // src/bcfy_grid.cpp / kk_grid.cpp, inputs/bcfy_eko/, inputs/kk_eko/ --
    // ported from inclusive-D0-UPC's QCDNUM-based bcfy_grid.cpp/kk_grid.cpp),
    // HymnD from the external prompt-D0 set.
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
