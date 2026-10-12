#include "amplitudelib.hpp"
#include "def.hpp"
#include <cmath>

using namespace std;

// Range of the logarithmic sampling of the gluon momentum k, up to pD0/z_h_min.
static void set_k_range(parameters* par)
{
    par->k_upper_frag = par->pD0 / par->z_h_min;
    par->log_k_frag = log(par->k_upper_frag / K_MIN);
}

double diffractiveCrossSection(void* p)
{
    parameters* par = (parameters*)p;
    set_k_range(par);

    gsl_monte_function F;
    F.f      = &integrand_diffractive;
    F.dim    = 7;
    F.params = par;

    // integration limits for {r1, r2, u_k, u_qp, b, u_xpo, z_h}
    bool eff = is_effective_flux(par);
    double b_lo = b_integration_min(par);
    double b_hi = eff ? 1.0 : par->bmax;
    double rmax = 99.0;
    double low[] = {0.,   0.,   0.,   0., b_lo, 0.,  par->z_h_min};
    double up[]  = {rmax, rmax, 1.0,  1., b_hi, 1.,  par->z_h_max};

    return vegas_integrate(&F, low, up, par->calls_diff);
}

double diffractiveCrossSection_xpom(void* p)
{
    parameters* par = (parameters*)p;
    set_k_range(par);

    gsl_monte_function F;
    F.f      = &integrand_diffractive_xpom;
    F.dim    = 6;
    F.params = par;

    // integration limits for {r1, r2, u_k, u_qp, b, z_h}
    bool eff = is_effective_flux(par);
    double b_lo = b_integration_min(par);
    double b_hi = eff ? 1.0 : par->bmax;
    double rmax = 99.0;
    double low[] = {0.,   0.,   0.,   0., b_lo, par->z_h_min};
    double up[]  = {rmax, rmax, 1.0,  1., b_hi, par->z_h_max};

    return vegas_integrate(&F, low, up, par->calls_diff);
}

// VEGAS integral of a fixed-q+ diffractive integrand: r1, r2 in [0, 99], the other variables in [0, 1].
static double vegas_fixed_qp(double (*integrand)(double*, size_t, void*), size_t dim, parameters* par)
{
    par->log_k = log(par->k_upper / K_MIN);

    gsl_monte_function F;
    F.f      = integrand;
    F.dim    = dim;
    F.params = par;

    double rmax = 99.0;
    double low[] = {0.,   0.,   0., 0.};
    double up[]  = {rmax, rmax, 1., 1.};

    return vegas_integrate(&F, low, up, par->calls_diff);
}

// Diffractive cross section at fixed q+, and at fixed q+ and x_po. Set par->k_upper before calling.
double diffractiveCrossSection_fixed_qp(void* p)
{
    return vegas_fixed_qp(&integrand_diffractive_fixed_qp, 4, (parameters*)p);
}

double diffractiveCrossSection_fixed_qp_xpom(void* p)
{
    return vegas_fixed_qp(&integrand_diffractive_fixed_qp_xpom, 3, (parameters*)p);
}
