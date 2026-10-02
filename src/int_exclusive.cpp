#include "amplitudelib.hpp"
#include "def.hpp"
#include <algorithm>
#include <cmath>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

using namespace std;


// Exclusive cross section: 5D integral below par->excl_pt_threshold,
// 3D integral at or above it.
double exclusiveCrossSection(void* p)
{
    parameters* par = (parameters*)p;
    bool factorized = par->pD0 >= par->excl_pt_threshold;
    size_t calls;
    if (factorized) calls = par->calls_excl_factorized;
    else             calls = par->calls_excl;

    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);

    gsl_monte_function F;
    F.params = par;

    // effective flux: b is not used, its range is [0,1]
    bool eff = is_effective_flux(par);
    double b_lo = b_integration_min(par);
    double b_hi = eff ? 1.0 : par->bmax;

    double low[5], up[5];
    if (factorized) {
        F.f   = &integrand_exclusive_factorized;
        F.dim = 3;
        // integration limits for {u_qp, b, zh}
        low[0] = 0.;    up[0] = 1.;
        low[1] = b_lo;  up[1] = b_hi;
        low[2] = par->zmin;  up[2] = par->zmax;
    } else {
        F.f   = &integrand_exclusive_mc;
        F.dim = 5;
        // integration limits for {r1, r2, u_qp, b, zh}
        double rmax = 99.0;
        low[0] = 0.;    up[0] = rmax;
        low[1] = 0.;    up[1] = rmax;
        low[2] = 0.;    up[2] = 1.;
        low[3] = b_lo;  up[3] = b_hi;
        low[4] = par->zmin;  up[4] = par->zmax;
    }

    double res, err;
    gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(F.dim);

    // short first run so VEGAS can adapt its grid
    gsl_monte_vegas_integrate(&F, low, up, F.dim, calls/10, r, s, &res, &err);

    // repeat until chi-square is close to 1, at most 3 times
    int iter = 0;
    do {
        gsl_monte_vegas_integrate(&F, low, up, F.dim, calls, r, s, &res, &err);
        iter++;
    } while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.25 && iter < 3);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    return res;
}


// Exclusive cross section at fixed x_P (par->fixed_xpo): VEGAS over {b, zh}.
double exclusiveCrossSection_xpom(void* p)
{
    parameters* par = (parameters*)p;
    size_t calls = par->calls_excl_factorized;

    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);

    gsl_monte_function F;
    F.f      = &integrand_exclusive_xpom;
    F.dim    = 2;
    F.params = par;

    bool eff = is_effective_flux(par);   // see exclusiveCrossSection() above
    double low[] = {b_integration_min(par), par->zmin};
    double up[]  = {eff ? 1.0 : par->bmax, par->zmax};

    double res, err;
    gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(F.dim);

    gsl_monte_vegas_integrate(&F, low, up, F.dim, calls/10, r, s, &res, &err);

    int iter = 0;
    do {
        gsl_monte_vegas_integrate(&F, low, up, F.dim, calls, r, s, &res, &err);
        iter++;
    } while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.25 && iter < 3);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    return res;
}


// Used for the two 1D integrals over r in exclusiveCrossSection_fixed_qp.
struct RadialParams {
    double pc;
    double m;
    int nu;
    double x_dip;
    AmplitudeLib* dipole;
};

static double radial_integrand(double r, void* params)
{
    RadialParams* rp = (RadialParams*)params;
    return r * Jn(rp->nu, rp->pc * r) * Kn(rp->nu, rp->m * r) * rp->dipole->N(r, rp->x_dip);
}

// Exclusive cross section at fixed q+ (q+ = 2p+, no flux, no fragmentation):
// two 1D integrals over r, I0 and I1. No Monte Carlo.
double exclusiveCrossSection_fixed_qp(void* p)
{
    parameters* par = (parameters*)p;

    double z1 = 0.5;   // q+ = 2p+, so p+/q+ = 1/2
    double z2 = 1.0 - z1;

    double pp   = (par->mt / M_SQRT2) * exp(par->y);
    double qp   = 2.0 * pp;
    double w2   = M_SQRT2 * qp * par->ss;
    double Mqq2 = (par->p2 + par->m2) / (z1 * z2);
    double xP   = Mqq2 / w2;
    if (xP <= 0.0 || xP > 0.1) return 0.0;
    double xP_dip = min(xP, 0.01);

    double rmax = 99.0;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

    RadialParams rp0;
    rp0.pc = par->p;
    rp0.m = par->m;
    rp0.nu = 0;
    rp0.x_dip = xP_dip;
    rp0.dipole = par->dipole;

    RadialParams rp1 = rp0;
    rp1.nu = 1;

    gsl_function F0;
    F0.function = &radial_integrand;
    F0.params = &rp0;

    gsl_function F1;
    F1.function = &radial_integrand;
    F1.params = &rp1;

    double I0, I1, err0, err1;
    gsl_integration_qags(&F0, 0.0, rmax, 1e-10, 1e-6, 1000, w, &I0, &err0);
    gsl_integration_qags(&F1, 0.0, rmax, 1e-10, 1e-6, 1000, w, &I1, &err1);

    gsl_integration_workspace_free(w);

    double factor = z1*z1 + z2*z2;
    double g = factor * I1*I1 + I0*I0;

    return z1 * par->m2 * g;
}
