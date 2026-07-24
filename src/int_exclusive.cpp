#include "amplitudelib.hpp"
#include "def.hpp"
#include <algorithm>
#include <cmath>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>

using namespace std;


double exclusiveCrossSection(void* p)
{
    parameters* par = (parameters*)p;
    size_t calls = par->calls;

    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);

    gsl_monte_function F;
    F.f      = &integrand_exclusive;
    F.dim    = 5;
    F.params = par;

    // {r1, r2, u_qp, b, zh} -- zh's box is the literal [zmin,zmax] range not [0,1])
    double rmax = 99.0;
    double low[] = {0.,   0.,   0., par->bmin, par->zmin};
    double up[]  = {rmax, rmax, 1., par->bmax, par->zmax};

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
