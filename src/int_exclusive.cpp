#include "amplitudelib.hpp"
#include "def.hpp"
#include <algorithm>
#include <cmath>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

using namespace std;


// Computes the exclusive cross section with a Monte Carlo integral over the
// integrand defined in integrand.cpp. Below par->excl_pt_threshold, uses the
// original, cheap-per-call 5D integrand (r1, r2 handed to VEGAS along with
// u_qp, b, zh); at or above it, switches to the factorized 3D one (r1, r2
// solved analytically per point, see integrand_exclusive_factorized), which
// converges properly at high pt but costs much more per call -- not worth
// paying for at low/mid pt, where the 5D version already converges fine.
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

    double low[5], up[5];
    if (factorized) {
        F.f   = &integrand_exclusive_factorized;
        F.dim = 3;
        // integration box for {u_qp, b, zh}
        // (zh's box is the actual [zmin, zmax] range, not [0,1] like u_qp)
        // r1, r2 are not here -- they're solved analytically inside the
        // integrand (see exclusive_radial_integrals in integrand.cpp).
        low[0] = 0.;         up[0] = 1.;
        low[1] = par->bmin;  up[1] = par->bmax;
        low[2] = par->zmin;  up[2] = par->zmax;
    } else {
        F.f   = &integrand_exclusive_mc;
        F.dim = 5;
        // integration box for {r1, r2, u_qp, b, zh}
        double rmax = 99.0;
        low[0] = 0.;    up[0] = rmax;
        low[1] = 0.;    up[1] = rmax;
        low[2] = 0.;    up[2] = 1.;
        low[3] = par->bmin;  up[3] = par->bmax;
        low[4] = par->zmin;  up[4] = par->zmax;
    }

    double res, err;
    gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(F.dim);

    // warm-up run so VEGAS can adapt its grid before the real integration
    gsl_monte_vegas_integrate(&F, low, up, F.dim, calls/10, r, s, &res, &err);

    // keep integrating until it converges (chi-square close to 1), or give up after 3 tries
    int iter = 0;
    do {
        gsl_monte_vegas_integrate(&F, low, up, F.dim, calls, r, s, &res, &err);
        iter++;
    } while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.25 && iter < 3);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    return res;
}


// Same as exclusiveCrossSection above, but at one fixed x_P (par->fixed_xpo)
// instead of integrating over the photon energy q+ -- see
// integrand_exclusive_xpom for why x_P being fixed eliminates q+ entirely
// (unlike diffractiveCrossSection_xpom, which still integrates u_qp/q+
// separately, since x_po there IS an independent variable). r1, r2 are
// solved analytically (exclusive_radial_integrals), same trick as the
// factorized high-pt integrand, so only {b, zh} are left for VEGAS.
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

    double low[] = {par->bmin, par->zmin};
    double up[]  = {par->bmax, par->zmax};

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


// Small helper struct + function used only by exclusiveCrossSection_fixedW
// below, to do the two simple 1D integrals over dipole size.
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

// "Fixed W, no photon flux" exclusive cross section. Here q+ = 3p+ is fixed
// by hand (z1 = 1/3) and pc is just par->p (no fragmentation step). Because
// par->p doesn't change during this calculation, the usual double integral
// over the two dipole sizes r1 and r2 splits into two separate, identical
// 1D integrals (I0 and I1). So instead of a Monte Carlo integration, we can
// just compute these two integrals directly with GSL and combine them.
double exclusiveCrossSection_fixedW(void* p)
{
    parameters* par = (parameters*)p;

    double z1 = 1.0 / 3.0;   // q+ = 3p+, so p+/q+ = 1/3
    double z2 = 1.0 - z1;

    double pp   = (par->mt / M_SQRT2) * exp(par->y);
    double qp   = 3.0 * pp;
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
