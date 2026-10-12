#include "def.hpp"
#include <cmath>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>

double vegas_integrate(gsl_monte_function* F, double* low, double* up, size_t calls)
{
    gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
    gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(F->dim);
    double res, err;

    // short first run so VEGAS can adapt its grid
    gsl_monte_vegas_integrate(F, low, up, F->dim, calls/10, r, s, &res, &err);

    // repeat until chi-square is close to 1, at most 3 times
    int iter = 0;
    do {
        gsl_monte_vegas_integrate(F, low, up, F->dim, calls, r, s, &res, &err);
        iter++;
    } while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.25 && iter < 3);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    return res;
}
