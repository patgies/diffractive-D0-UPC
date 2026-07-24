#include "amplitudelib.hpp"
#include "def.hpp"
#include <cmath>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>

using namespace std;


double diffractiveCrossSection(void* p)
{
    parameters* par = (parameters*)p;
    size_t calls = par->calls;

    constexpr double k_lo_frag = 0.01;   
    par->k_upper_frag = par->pD0 / par->zmin;
    par->log_k_frag = log(par->k_upper_frag / k_lo_frag);

    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);

    gsl_monte_function F;
    F.f      = &integrand_diffractive;
    F.dim    = 7;
    F.params = par;

    // {r1, r2, u_k, u_qp, b, u_xpo, zh}; u_k, u_xpo in [0,1], log-mapped;
    // zh's box is the literal [zmin,zmax] range (not [0,1])
    double rmax = 99.0;
    double low[] = {0.,   0.,   0.,   0., par->bmin, 0.,  par->zmin};
    double up[]  = {rmax, rmax, 1.0,  1., par->bmax, 1.,  par->zmax};

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
