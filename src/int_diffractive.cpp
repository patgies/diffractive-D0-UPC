#include "amplitudelib.hpp"
#include "def.hpp"
#include "alphas_running.hpp"
#include <cmath>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>

using namespace std;


// Diffractive cross section: 7D VEGAS integral.
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

    // integration limits for {r1, r2, u_k, u_qp, b, u_xpo, zh}
    bool eff = is_effective_flux(par);
    double b_lo = b_integration_min(par);
    double b_hi = eff ? 1.0 : par->bmax;
    double rmax = 99.0;
    double low[] = {0.,   0.,   0.,   0., b_lo, 0.,  par->zmin};
    double up[]  = {rmax, rmax, 1.0,  1., b_hi, 1.,  par->zmax};

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

// Diffractive cross section at fixed x_po (par->fixed_xpo): 6D VEGAS integral.
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

    // integration limits for {r1, r2, u_k, u_qp, b, zh}
    bool eff = is_effective_flux(par);
    double b_lo = b_integration_min(par);
    double b_hi = eff ? 1.0 : par->bmax;
    double rmax = 99.0;
    double low[] = {0.,   0.,   0.,   0., b_lo, par->zmin};
    double up[]  = {rmax, rmax, 1.0,  1., b_hi, par->zmax};

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

// Diffractive cross section at fixed q+: 4D VEGAS integral over {r1, r2, u_k, u_xpo}.
// The result has the full prefactor. Set par->k_upper before calling.
double diffractiveCrossSection_fixed_qp(void* p)
{
    parameters* par = (parameters*)p;
    size_t calls = par->calls_diff;

    double k_lo = 0.01;
    par->log_k = log(par->k_upper / k_lo);

    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);

    gsl_monte_function F;
    F.f      = &integrand_diffractive_fixed_qp;
    F.dim    = 4;
    F.params = par;

    // integration limits for {r1, r2, u_k, u_xpo}
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

    // prefactor: alpha_s * alpha_em * e_c^2 * (Nc^2-1) * sigma0 / (8*pi^4)
    double alphae = 1.0 / 137.0;
    double e_c    = 2.0 / 3.0;
    double Nc     = 3.0;
    double sigma0 = 16.36;   // dipole normalization, same value as in the python scripts
    double prefactor = alphas_run(par->mt) * alphae * e_c*e_c * (Nc*Nc - 1.0) * sigma0
                        / (8.0 * M_PI*M_PI*M_PI*M_PI);

    return res * prefactor;
}
