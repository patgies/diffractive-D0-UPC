#include "amplitudelib.hpp"
#include "def.hpp"
#include <algorithm>
#include <cmath>

using namespace std;


double exclusiveCrossSection(void* p)
{
    parameters* par = (parameters*)p;
    bool factorized = par->pD0 >= par->excl_pt_threshold;

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
        // integration limits for {u_qp, b, z_h}
        low[0] = 0.;    up[0] = 1.;
        low[1] = b_lo;  up[1] = b_hi;
        low[2] = par->z_h_min;  up[2] = par->z_h_max;
    } else {
        F.f   = &integrand_exclusive_mc;
        F.dim = 5;
        // integration limits for {r1, r2, u_qp, b, z_h}
        double rmax = 99.0;
        low[0] = 0.;    up[0] = rmax;
        low[1] = 0.;    up[1] = rmax;
        low[2] = 0.;    up[2] = 1.;
        low[3] = b_lo;  up[3] = b_hi;
        low[4] = par->z_h_min;  up[4] = par->z_h_max;
    }

    return vegas_integrate(&F, low, up, factorized ? par->calls_excl_factorized : par->calls_excl);
}


double exclusiveCrossSection_xpom(void* p)
{
    parameters* par = (parameters*)p;

    gsl_monte_function F;
    F.f      = &integrand_exclusive_xpom;
    F.dim    = 2;
    F.params = par;

    bool eff = is_effective_flux(par);
    double low[] = {b_integration_min(par), par->z_h_min};
    double up[]  = {eff ? 1.0 : par->bmax, par->z_h_max};

    return vegas_integrate(&F, low, up, par->calls_excl_factorized);
}


// Exclusive cross section at fixed q+ (q+ = 2p+, no flux, no fragmentation):
// two 1D integrals over r, I0 and I1. No Monte Carlo.
double exclusiveCrossSection_fixed_qp(void* p)
{
    parameters* par = (parameters*)p;

    double z = 0.5;   // q+ = 2p+, so p+/q+ = 1/2
    double zbar = 1.0 - z;

    double pp   = (par->mt / M_SQRT2) * exp(par->y);
    double qp   = 2.0 * pp;
    double w2   = M_SQRT2 * qp * par->ss;
    double Mqq2 = (par->p2 + par->m2) / (z * zbar);
    double xP   = Mqq2 / w2;
    if (xP <= 0.0 || xP > XPO_MAX) return 0.0;
    double xP_dip = min(xP, XPO_FROZEN);

    double I0, I1;
    exclusive_radial_integrals(par->p, par->m, xP_dip, par->dipole, I0, I1);

    double factor = z*z + zbar*zbar;
    double g = factor * I1*I1 + I0*I0;

    return z * par->m2 * g;
}
