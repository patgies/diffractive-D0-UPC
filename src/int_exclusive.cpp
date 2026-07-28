#include "amplitudelib.hpp"
#include "def.hpp"
#include <algorithm>
#include <cmath>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

using namespace std;


// Computes the exclusive cross section by doing a 5-dimensional Monte Carlo
// integral over the integrand defined in integrand.cpp.
double exclusiveCrossSection(void* p)
{
    parameters* par = (parameters*)p;
    size_t calls = par->calls_excl;

    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);

    gsl_monte_function F;
    F.f      = &integrand_exclusive;
    F.dim    = 5;
    F.params = par;

    // integration box for {r1, r2, u_qp, b, zh}
    // (zh's box is the actual [zmin, zmax] range, not [0,1] like the others)
    double rmax = 99.0;
    double low[] = {0.,   0.,   0., par->bmin, par->zmin};
    double up[]  = {rmax, rmax, 1., par->bmax, par->zmax};

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
