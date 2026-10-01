#include "def.hpp"
#include "amplitudelib.hpp"
#include "fragmentation.hpp"
#include <cmath>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

using namespace std;

double Kn(int nu, double x)
{
    gsl_sf_result result;
    int status = gsl_sf_bessel_Kn_e(nu, x, &result);
    if (status != GSL_SUCCESS) return 0.0;
    return result.val;
}

// Bessel function of the first kind, J_nu(x). Same idea as Kn above.
double Jn(int nu, double x)
{
    gsl_sf_result result;
    int status = gsl_sf_bessel_Jn_e(nu, x, &result);
    if (status != GSL_SUCCESS) return 0.0;
    return result.val;
}


// Fragmentation function D(zh)/zh^2, for turning a charm quark into a D0 meson.
// Picks the right formula depending on which fragmentation model we're using.
static double fragmentation(double zh, parameters* par)
{
    // All three frag_types (BCFY, KniehlKramer, HymnD) are now DGLAP-evolved
    // z-interpolators built once per point at frag_scale=Q (see main.cpp) --
    // BCFY/KniehlKramer from their own eko-evolved grids (src/bcfy_grid.cpp,
    // src/kk_grid.cpp), HymnD from the external prompt-D0 set. The raw,
    // non-evolved Dc_to_D0()/D_kniehl_kramer() (fragmentation.cpp) are still
    // used as the mu=mc seed for that evolution, just not called directly
    // here anymore.
    double D_frag = par->D_frag_interp->Evaluate(zh);
    return D_frag / (zh * zh);
}

// Radial (r) integrand for the exclusive process: r * J_nu(pc*r) * K_nu(m*r) * N(r, x_dip),
// for one dipole size. Used by exclusive_radial_integrals below.
struct ExclRadialParams {
    double pc, m, x_dip;
    int nu;
    AmplitudeLib* dipole;
};

static double excl_radial_integrand(double r, void* params)
{
    ExclRadialParams* rp = (ExclRadialParams*)params;
    return r * Jn(rp->nu, rp->pc * r) * Kn(rp->nu, rp->m * r) * rp->dipole->N(r, rp->x_dip);
}

// For fixed pc, m and x_dip, the exclusive integrand's r1 and r2 dependence
// factorizes into two separate 1D integrals I0, I1 (same trick as
// exclusiveCrossSection_fixed_qp in int_exclusive.cpp). Solving them here with
// deterministic adaptive quadrature avoids handing r1, r2 to the Monte Carlo
// routine, where the oscillating Bessel functions converge very slowly at
// high pt. 
static void exclusive_radial_integrals(double pc, double m, double x_dip, AmplitudeLib* dipole, double& I0, double& I1)
{
    static thread_local gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    const double rmax = 99.0;

    ExclRadialParams rp0{pc, m, x_dip, 0, dipole};
    ExclRadialParams rp1{pc, m, x_dip, 1, dipole};

    gsl_function F0; F0.function = &excl_radial_integrand; F0.params = &rp0;
    gsl_function F1; F1.function = &excl_radial_integrand; F1.params = &rp1;

    double err0, err1;
    gsl_integration_qags(&F0, 0.0, rmax, 1e-10, 1e-6, 1000, w, &I0, &err0);
    gsl_integration_qags(&F1, 0.0, rmax, 1e-10, 1e-6, 1000, w, &I1, &err1);
}

// Exclusive integrand at the D0 level: dsigma/(d2pD0 dy), factorized version.
// This function gets called many times by the Monte Carlo integration
// routine, once for each random point "vec" it picks. r1 and r2 are not
// among the 3 MC variables -- they're solved analytically inside via
// exclusive_radial_integrals; only {u_qp, b, zh} are left for
// VEGAS. 
double integrand_exclusive_factorized(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    // the 3 numbers the integration routine gives us this time
    double u_qp = vec[0];
    double b    = vec[1];
    double zh   = vec[2];

    double pc = par->pD0 / zh;   // charm quark's transverse momentum
    double m  = par->m;

    double p2 = pc * pc;
    double mt = sqrt(p2 + par->m2);
    double pp = (mt / M_SQRT2) * exp(par->y);

    // turn u_qp (between 0 and 1) into qp (between pp and qpmax)
    double jac_qp = par->qpmax - pp;
    if (jac_qp <= 0.0) return 0.0;
    double qp = pp + u_qp * jac_qp;

    // par->flux_model=="EFF": the whole b-dependence (and target-nucleus
    // geometric convolution, arXiv:2404.09731 Eq. 4) is already integrated
    // out into effective_photon_flux(qp,p); the sampled "b" here is then an
    // unused placeholder (its VEGAS box is set to [0,1] in int_exclusive.cpp/
    // int_diffractive.cpp for this mode, contributing a no-op Jacobian).
    double flux = is_effective_flux(par) ? effective_photon_flux(qp, p) : photon_flux(b, qp, p);

    // light-cone momentum fractions
    double z1 = pp / qp;
    double z2 = 1.0 - z1;

    // photon-proton energy squared (W^2), quark-pair mass squared, and
    // the momentum fraction of the pomeron
    double w2   = M_SQRT2 * qp * par->ss;
    double Mqq2 = (p2 + par->m2) / (z1 * z2);
    double xP   = Mqq2 / w2;
    if (xP <= 0.0 || xP > 0.1) return 0.0;
    double xP_dip = min(xP, 0.01);

    double I0, I1;
    exclusive_radial_integrals(pc, m, xP_dip, par->dipole, I0, I1);

    double factor = z1*z1 + z2*z2;
    double g = factor * I1*I1 + I0*I0;

    double kernel = z1 * m*m * jac_qp * flux * g;

    return kernel * fragmentation(zh, par);
}

// Exclusive integrand at the D0 level, at one FIXED x_P (par->fixed_xpo)
// instead of integrating over q+ (u_qp dropped). Unlike the diffractive
// channel's x_po, x_P is not an independent variable here -- it's a function
// of q+ (and y, mt) alone, so "fixing x_P" means picking out the single q+
// that gives it (a delta function, not holding one of several independent
// variables fixed), and multiplying by the resulting Jacobian.
//
// Solve for z_* = p+/q+_* and q+_* at fixed x_P, y (mt depends on zh, so
// this is solved fresh per zh sample):
//   x_P = m_t/(sqrt(s_NN)) * e^-y / (1-z_*)     [see notes/ for derivation;
//                                                 corrects a lost reciprocal
//                                                 relative to the pasted
//                                                 source formula]
//   q+_* = p+ / z_*
// Then (see e.g. Eq. 25 of the reference this was ported from):
//   dsigma/(dy d2p dlnx_P) = dsigma/(dy d2p dq+)|_{q+=q+_*} * (q+_*)^2(1-z_*)/p+
// To keep the SAME raw-output convention as integrand_diffractive_xpom
// (plain x_P density, dsigma/dx_P -- see xpom_spectrum.py's "xpo *
// dsigma_dxpo" ln-conversion step, applied uniformly to both channels),
// this returns dsigma/dx_P, i.e. that Jacobian divided by one more power
// of x_P.
double integrand_exclusive_xpom(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double b  = vec[0];
    double zh = vec[1];

    double xP = par->fixed_xpo;
    if (xP <= 0.0 || xP > 0.1) return 0.0;

    double pc = par->pD0 / zh;
    double m  = par->m;
    double p2 = pc * pc;
    double mt = sqrt(p2 + par->m2);
    double pp = (mt / M_SQRT2) * exp(par->y);   // p+

    double z2 = mt * exp(-par->y) / (par->ss * xP);   // 1 - z_*
    double z1 = 1.0 - z2;
    if (z1 <= 0.0 || z2 <= 0.0) return 0.0;   // this (pD0, y, zh) can't reach this x_P

    double qp = pp / z1;   // q+_*, Eq. 19
    if (qp > par->qpmax) return 0.0;

    // par->flux_model=="EFF": the whole b-dependence (and target-nucleus
    // geometric convolution, arXiv:2404.09731 Eq. 4) is already integrated
    // out into effective_photon_flux(qp,p); the sampled "b" here is then an
    // unused placeholder (its VEGAS box is set to [0,1] in int_exclusive.cpp/
    // int_diffractive.cpp for this mode, contributing a no-op Jacobian).
    double flux = is_effective_flux(par) ? effective_photon_flux(qp, p) : photon_flux(b, qp, p);

    double xP_dip = min(xP, 0.01);
    double I0, I1;
    exclusive_radial_integrals(pc, m, xP_dip, par->dipole, I0, I1);

    double factor = z1*z1 + z2*z2;
    double g = factor * I1*I1 + I0*I0;

    // (q+_*)^2 (1-z_*) / p+, divided by x_P (see docstring above)
    double jac_xpo = qp*qp * z2 / (pp * xP);

    double kernel = z1 * m*m * jac_xpo * flux * g;

    return kernel * fragmentation(zh, par);
}

// Exclusive integrand at the D0 level: dsigma/(d2pD0 dy), original 5D
// version (r1, r2 handed to VEGAS along with u_qp, b, zh instead of being
// solved analytically). At high pt, pc = pD0/zh gets large and Jn(pc*r)
// oscillates fast, which VEGAS's per-axis grid adaptation cannot resolve.  Kept
// around because it converges fine at low/mid pt,
// where exclusiveCrossSection still uses it.
double integrand_exclusive_mc(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double r1   = vec[0];
    double r2   = vec[1];
    double u_qp = vec[2];
    double b    = vec[3];
    double zh   = vec[4];

    double pc = par->pD0 / zh;   // charm quark's transverse momentum
    double m  = par->m;

    double p2 = pc * pc;
    double mt = sqrt(p2 + par->m2);
    double pp = (mt / M_SQRT2) * exp(par->y);

    // turn u_qp (between 0 and 1) into qp (between pp and qpmax)
    double jac_qp = par->qpmax - pp;
    if (jac_qp <= 0.0) return 0.0;
    double qp = pp + u_qp * jac_qp;

    // par->flux_model=="EFF": the whole b-dependence (and target-nucleus
    // geometric convolution, arXiv:2404.09731 Eq. 4) is already integrated
    // out into effective_photon_flux(qp,p); the sampled "b" here is then an
    // unused placeholder (its VEGAS box is set to [0,1] in int_exclusive.cpp/
    // int_diffractive.cpp for this mode, contributing a no-op Jacobian).
    double flux = is_effective_flux(par) ? effective_photon_flux(qp, p) : photon_flux(b, qp, p);

    // light-cone momentum fractions
    double z1 = pp / qp;
    double z2 = 1.0 - z1;

    // photon-proton energy squared (W^2), quark-pair mass squared, and
    // the momentum fraction of the pomeron
    double w2   = M_SQRT2 * qp * par->ss;
    double Mqq2 = (p2 + par->m2) / (z1 * z2);
    double xP   = Mqq2 / w2;
    if (xP <= 0.0 || xP > 0.1) return 0.0;
    double xP_dip = min(xP, 0.01);

    // dipole scattering amplitude for each of the two dipole sizes
    double D1 = par->dipole->N(r1, xP_dip);
    double D2 = par->dipole->N(r2, xP_dip);

    double f1_r1 = r1 * Jn(1, pc*r1) * Kn(1, m*r1) * D1;
    double f1_r2 = r2 * Jn(1, pc*r2) * Kn(1, m*r2) * D2;
    double f0_r1 = r1 * Jn(0, pc*r1) * Kn(0, m*r1) * D1;
    double f0_r2 = r2 * Jn(0, pc*r2) * Kn(0, m*r2) * D2;

    double factor = z1*z1 + z2*z2;
    double g = factor * f1_r1 * f1_r2 + f0_r1 * f0_r2;

    double kernel = z1 * m*m * jac_qp * flux * g;

    return kernel * fragmentation(zh, par);
}

// Diffractive integrand at the D0 level: dsigma/(d2pD0 dy), summed over the
// pomeron momentum fraction x_po and the emitted gluon's momentum k.
double integrand_diffractive(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double k_lo_frag = 0.01;
    double x_po_min_log = 1e-6;

    // the 7 numbers the integration routine gives us this time
    double r1   = vec[0];
    double r2   = vec[1];
    double k    = k_lo_frag * exp(vec[2] * par->log_k_frag);
    double u_qp = vec[3];
    double b    = vec[4];
    double log_xpo_range = log(0.1 / x_po_min_log);
    double x_po = x_po_min_log * exp(vec[5] * log_xpo_range);
    double zh   = vec[6];

    double pc = par->pD0 / zh;
    if (k > pc) return 0.0;   // the gluon can't be more energetic than the charm quark that radiated it

    double m2 = par->m2;
    double p2 = pc * pc;
    double mt = sqrt(p2 + m2);
    double pp = (mt / M_SQRT2) * exp(par->y);

    // turn u_qp (between 0 and 1) into qp (between pp and qpmax)
    double jac_qp = par->qpmax - pp;
    if (jac_qp <= 0.0) return 0.0;
    double qp = pp + u_qp * jac_qp;

    // par->flux_model=="EFF": the whole b-dependence (and target-nucleus
    // geometric convolution, arXiv:2404.09731 Eq. 4) is already integrated
    // out into effective_photon_flux(qp,p); the sampled "b" here is then an
    // unused placeholder (its VEGAS box is set to [0,1] in int_exclusive.cpp/
    // int_diffractive.cpp for this mode, contributing a no-op Jacobian).
    double flux = is_effective_flux(par) ? effective_photon_flux(qp, p) : photon_flux(b, qp, p);

    // light-cone momentum fractions
    double z1 = pp / qp;
    double z2 = 1.0 - z1;

    // photon-proton energy squared (W^2), quark-pair mass squared, 
    double w2   = M_SQRT2 * qp * par->ss;
    double Mqq2 = (p2 + m2) / (z1 * z2);
    double x    = Mqq2 / (x_po * w2);
    if (x <= 0.0 || x >= 1.0) return 0.0;

    double xfac    = x / (1.0 - x);
    double omega2  = xfac * k * k;
    double sqomega = sqrt(xfac) * k;

    // the "hard factor" H^2
    double denom4      = (p2 + m2) * (p2 + m2);
    double p4m4        = p2*p2 + m2*m2;
    double two_m2p2m2  = 2.0 * m2 * p2;
    double H2 = z1 / (denom4 * denom4) * ((z1*z1 + z2*z2) * p4m4 + two_m2p2m2);

    // the gluon transverse-momentum-dependent gDTMD^2
    double gDTMD2 = omega2 * omega2 * r1*r2
                    * Jn(2, k*r1) * Kn(2, sqomega*r1)
                    * Jn(2, k*r2) * Kn(2, sqomega*r2);

    double fun = H2 * (1.0 / (1.0 - x)) * gDTMD2;

    // dipole scattering amplitude for each of the two dipole sizes
    double x_po_dip = min(x_po, 0.01);
    double D1 = par->dipole->N(r1, x_po_dip);
    double D2 = par->dipole->N(r2, x_po_dip);

    double kernel = 2.0*M_PI * k * jac_qp * flux * D1*D2 * fun
                    * (k * par->log_k_frag) * (x_po * log_xpo_range);

    return kernel * fragmentation(zh, par);
}

// Same physics as integrand_diffractive above, but x_po is fixed at
// par->fixed_xpo instead of being one of the random numbers the integration
// routine picks. So this gives dsigma/(d2pD0 dy dx_po) at one specific
// x_po, integrated over the remaining 6 variables {r1, r2, u_k, u_qp, b, zh}.
double integrand_diffractive_xpom(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double k_lo_frag = 0.01;

    double r1   = vec[0];
    double r2   = vec[1];
    double k    = k_lo_frag * exp(vec[2] * par->log_k_frag);
    double u_qp = vec[3];
    double b    = vec[4];
    double zh   = vec[5];
    double x_po = par->fixed_xpo;

    double pc = par->pD0 / zh;
    if (k > pc) return 0.0;   

    double m2 = par->m2;
    double p2 = pc * pc;
    double mt = sqrt(p2 + m2);
    double pp = (mt / M_SQRT2) * exp(par->y);

    // turn u_qp (between 0 and 1) into qp (between pp and qpmax)
    double jac_qp = par->qpmax - pp;
    if (jac_qp <= 0.0) return 0.0;
    double qp = pp + u_qp * jac_qp;

    // par->flux_model=="EFF": the whole b-dependence (and target-nucleus
    // geometric convolution, arXiv:2404.09731 Eq. 4) is already integrated
    // out into effective_photon_flux(qp,p); the sampled "b" here is then an
    // unused placeholder (its VEGAS box is set to [0,1] in int_exclusive.cpp/
    // int_diffractive.cpp for this mode, contributing a no-op Jacobian).
    double flux = is_effective_flux(par) ? effective_photon_flux(qp, p) : photon_flux(b, qp, p);

    // light-cone momentum fractions
    double z1 = pp / qp;
    double z2 = 1.0 - z1;

    // photon-proton energy squared (W^2), quark-pair mass squared,
    double w2   = M_SQRT2 * qp * par->ss;
    double Mqq2 = (p2 + m2) / (z1 * z2);
    double x    = Mqq2 / (x_po * w2);
    if (x <= 0.0 || x >= 1.0) return 0.0;

    double xfac    = x / (1.0 - x);
    double omega2  = xfac * k * k;
    double sqomega = sqrt(xfac) * k;

    // the "hard factor" H^2
    double denom4      = (p2 + m2) * (p2 + m2);
    double p4m4        = p2*p2 + m2*m2;
    double two_m2p2m2  = 2.0 * m2 * p2;
    double H2 = z1 / (denom4 * denom4) * ((z1*z1 + z2*z2) * p4m4 + two_m2p2m2);

    // the gluon transverse-momentum-dependent gDTMD^2
    double gDTMD2 = omega2 * omega2 * r1*r2
                    * Jn(2, k*r1) * Kn(2, sqomega*r1)
                    * Jn(2, k*r2) * Kn(2, sqomega*r2);

    double fun = H2 * (1.0 / (1.0 - x)) * gDTMD2;

    // dipole scattering amplitude for each of the two dipole sizes
    double x_po_dip = min(x_po, 0.01);
    double D1 = par->dipole->N(r1, x_po_dip);
    double D2 = par->dipole->N(r2, x_po_dip);

    double kernel = 2.0*M_PI * k * jac_qp * flux * D1*D2 * fun
                    * (k * par->log_k_frag);

    return kernel * fragmentation(zh, par);
}

// "Fixed q+, no photon flux" diffractive integrand. Here q+ is set to 2*p+
// by hand (so z1 = 1/2 always) and pc is just par->p directly. There's no
// fragmentation step. Only 4 random numbers this time: {r1, r2, u_k, u_xpo}.

double integrand_diffractive_fixed_qp(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double k_lo = 0.01;
    double x_po_min_log = 1e-6;

    double r1   = vec[0];
    double r2   = vec[1];
    double k    = k_lo * exp(vec[2] * par->log_k);
    double log_xpo_range = log(0.1 / x_po_min_log);
    double x_po = x_po_min_log * exp(vec[3] * log_xpo_range);

    if (k > par->p) return 0.0;   

    double z1 = 0.5;   // q+ = 2p+, so p+/q+ = 1/2
    double z2 = 1.0 - z1;

    double pp = (par->mt / M_SQRT2) * exp(par->y);
    double qp = 2.0 * pp;
    double w2   = M_SQRT2 * qp * par->ss;
    double Mqq2 = (par->p2 + par->m2) / (z1 * z2);
    double x    = Mqq2 / (x_po * w2);
    if (x <= 0.0 || x >= 1.0) return 0.0;

    double xfac    = x / (1.0 - x);
    double omega2  = xfac * k * k;
    double sqomega = sqrt(xfac) * k;

    // the "hard factor" H^2,
    double H2 = z1 * par->inv_denom8 * ((z1*z1 + z2*z2) * par->p4m4 + par->two_m2p2m2);

    // the gluon transverse-momentum-dependent gDTMD^2
    double gDTMD2 = omega2 * omega2 * r1*r2
                    * Jn(2, k*r1) * Kn(2, sqomega*r1)
                    * Jn(2, k*r2) * Kn(2, sqomega*r2);

    double fun = H2 * (1.0 / (1.0 - x)) * gDTMD2;

    // dipole scattering amplitude for each of the two dipole sizes
    double x_po_dip = min(x_po, 0.01);
    double D1 = par->dipole->N(r1, x_po_dip);
    double D2 = par->dipole->N(r2, x_po_dip);

    double kernel = 2.0*M_PI * k * D1*D2 * fun
                    * (k * par->log_k) * (x_po * log_xpo_range);

    return kernel;
}

// Same as integrand_diffractive_fixed_qp above, but x_po is fixed at
// par->fixed_xpo instead of being one of the random integration
// variables (dropping u_xpo and its log-mapping Jacobian) -- gives
// dsigma_fixed_qp/(d2K dx_po) at one specific x_po, integrated over the
// remaining 3 variables {r1, r2, u_k}.
double integrand_diffractive_fixed_qp_xpom(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double k_lo = 0.01;

    double r1   = vec[0];
    double r2   = vec[1];
    double k    = k_lo * exp(vec[2] * par->log_k);
    double x_po = par->fixed_xpo;

    if (k > par->p) return 0.0;

    double z1 = 1.0 / 3.0;   // q+ = 3p+, so p+/q+ = 1/3
    double z2 = 1.0 - z1;

    double pp = (par->mt / M_SQRT2) * exp(par->y);
    double qp = 3.0 * pp;
    double w2   = M_SQRT2 * qp * par->ss;
    double Mqq2 = (par->p2 + par->m2) / (z1 * z2);
    double x    = Mqq2 / (x_po * w2);
    if (x <= 0.0 || x >= 1.0) return 0.0;

    double xfac    = x / (1.0 - x);
    double omega2  = xfac * k * k;
    double sqomega = sqrt(xfac) * k;

    // the "hard factor" H^2,
    double H2 = z1 * par->inv_denom8 * ((z1*z1 + z2*z2) * par->p4m4 + par->two_m2p2m2);

    // the gluon transverse-momentum-dependent gDTMD^2
    double gDTMD2 = omega2 * omega2 * r1*r2
                    * Jn(2, k*r1) * Kn(2, sqomega*r1)
                    * Jn(2, k*r2) * Kn(2, sqomega*r2);

    double fun = H2 * (1.0 / (1.0 - x)) * gDTMD2;

    // dipole scattering amplitude for each of the two dipole sizes
    double x_po_dip = min(x_po, 0.01);
    double D1 = par->dipole->N(r1, x_po_dip);
    double D2 = par->dipole->N(r2, x_po_dip);

    double kernel = 2.0*M_PI * k * D1*D2 * fun
                    * (k * par->log_k);

    return kernel;
}
