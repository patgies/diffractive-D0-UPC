#include "amplitudelib.hpp"
#include "def.hpp"
#include "alphas_running.hpp"
#include <cmath>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>

using namespace std;


// Computes the diffractive cross section with a 7-dimensional Monte Carlo
// integral over the integrand defined in integrand.cpp.
double diffractiveCrossSection(void* p)
{
    parameters* par = (parameters*)p;
    size_t calls = par->calls_diff;

    double k_lo_frag = 0.01;
    par->k_upper_frag = par->pD0 / par->zmin;
    par->log_k_frag = log(par->k_upper_frag / k_lo_frag);

    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);

    gsl_monte_function F;
    F.f      = &integrand_diffractive;
    F.dim    = 7;
    F.params = par;

    // integration box for {r1, r2, u_k, u_qp, b, u_xpo, zh}
    // (u_k and u_xpo are between 0 and 1 and get log-mapped inside the
    // integrand; zh's box is the actual [zmin, zmax] range)
    double rmax = 99.0;
    double low[] = {0.,   0.,   0.,   0., par->bmin, 0.,  par->zmin};
    double up[]  = {rmax, rmax, 1.0,  1., par->bmax, 1.,  par->zmax};

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

// Same as diffractiveCrossSection above, but x_po is fixed at par->fixed_xpo
// instead of being one of the random integration variables. So this is a
// 6-dimensional integral over {r1, r2, u_k, u_qp, b, zh}, and gives
// dsigma/(d2pD0 dy dx_po) at that one value of x_po.
double diffractiveCrossSection_xpom(void* p)
{
    parameters* par = (parameters*)p;
    size_t calls = par->calls_diff;

    double k_lo_frag = 0.01;
    par->k_upper_frag = par->pD0 / par->zmin;
    par->log_k_frag = log(par->k_upper_frag / k_lo_frag);

    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);

    gsl_monte_function F;
    F.f      = &integrand_diffractive_xpom;
    F.dim    = 6;
    F.params = par;

    // integration box for {r1, r2, u_k, u_qp, b, zh}
    double rmax = 99.0;
    double low[] = {0.,   0.,   0.,   0., par->bmin, par->zmin};
    double up[]  = {rmax, rmax, 1.0,  1., par->bmax, par->zmax};

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

// "Fixed W, no photon flux" diffractive cross section: a 4-dimensional
// integral over {r1, r2, u_k, u_xpo}.
// There's no {u_qp, b, zh} here because q+ = 3p+ is fixed by hand, there's
// no photon flux to integrate over, and pc is just par->p directly (no
// fragmentation). par->k_upper needs to already be set by the caller.
//
// Unlike diffractiveCrossSection/diffractiveCrossSection_xpom above, this
// one multiplies in the full physical prefactor at the end (including
// alpha_s), that's possible here because mT (used for alpha_s) is fixed for this whole
// calculation, rather than changing every time zh is sampled.
double diffractiveCrossSection_fixedW(void* p)
{
    parameters* par = (parameters*)p;
    size_t calls = par->calls_diff;

    double k_lo = 0.01;
    par->log_k = log(par->k_upper / k_lo);

    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);

    gsl_monte_function F;
    F.f      = &integrand_diffractive_fixedW;
    F.dim    = 4;
    F.params = par;

    // integration box for {r1, r2, u_k, u_xpo}
    double rmax = 99.0;
    double low[] = {0.,   0.,   0., 0.};
    double up[]  = {rmax, rmax, 1., 1.};

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

    // physics prefactor: alpha_s * alpha_em * e_c^2 * (Nc^2-1) * sigma0 / (8*pi^4)
    double alphae = 1.0 / 137.0;
    double e_c    = 2.0 / 3.0;
    double Nc     = 3.0;
    double sigma0 = 16.36;   // dipole normalization, same value used in the python scripts
    double prefactor = alphas_run(par->mt) * alphae * e_c*e_c * (Nc*Nc - 1.0) * sigma0
                        / (8.0 * M_PI*M_PI*M_PI*M_PI);

    return res * prefactor;
}
