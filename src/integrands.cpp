#include "def.hpp"
#include "amplitudelib.hpp"
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

double Jn(int nu, double x)
{
    gsl_sf_result result;
    int status = gsl_sf_bessel_Jn_e(nu, x, &result);
    if (status != GSL_SUCCESS) return 0.0;
    return result.val;
}


// Fragmentation factor D(z_h)/z_h^2 for c -> D0.
static double fragmentation(double z_h, parameters* par)
{
    double D_frag = par->D_frag_interp->Evaluate(z_h);
    return D_frag / (z_h * z_h);
}

// Function of r for the exclusive process: r * J_nu(pc*r) * K_nu(m*r) * N(r, x_dip).
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

// At fixed pc, m and x_dip the r1 and r2 integrals split into two 1D integrals, I0 and I1,
// computed here with GSL.
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

// Exclusive, 3D version: VEGAS integral over {u_qp, b, z_h}. r1 and r2 are done inside.
double integrand_exclusive_factorized(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double u_qp = vec[0];
    double b    = vec[1];
    double z_h   = vec[2];

    double pc = par->pD0 / z_h;   // charm quark's transverse momentum
    double m  = par->m;

    double p2 = pc * pc;
    double mt = sqrt(p2 + par->m2);
    double pp = (mt / M_SQRT2) * exp(par->y);

    // change u_qp (0 to 1) into qp (pp to qpmax)
    double jac_qp = par->qpmax - pp;
    if (jac_qp <= 0.0) return 0.0;
    double qp = pp + u_qp * jac_qp;

    // effective flux: b is not used
    double flux = is_effective_flux(par) ? effective_photon_flux(qp, p) : photon_flux(b, qp, p);

    // light-cone momentum fractions
    double z = pp / qp;
    double zbar = 1.0 - z;

    // W^2, quark-pair mass squared and pomeron momentum fraction
    double w2   = M_SQRT2 * qp * par->ss;
    double Mqq2 = (p2 + par->m2) / (z * zbar);
    double xP   = Mqq2 / w2;
    if (xP <= 0.0 || xP > 0.1) return 0.0;
    double xP_dip = min(xP, 0.01);

    double I0, I1;
    exclusive_radial_integrals(pc, m, xP_dip, par->dipole, I0, I1);

    double factor = z*z + zbar*zbar;
    double g = factor * I1*I1 + I0*I0;

    double kernel = z * m*m * jac_qp * flux * g;

    return kernel * fragmentation(z_h, par);
}

// Exclusive at fixed x_P (par->fixed_xpo): x_P fixes q+, so VEGAS only
// integrates over {b, z_h}. Returns dsigma/dx_P.
double integrand_exclusive_xpom(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double b  = vec[0];
    double z_h = vec[1];

    double xP = par->fixed_xpo;
    if (xP <= 0.0 || xP > 0.1) return 0.0;

    double pc = par->pD0 / z_h;
    double m  = par->m;
    double p2 = pc * pc;
    double mt = sqrt(p2 + par->m2);
    double pp = (mt / M_SQRT2) * exp(par->y);   // p+

    double zbar = mt * exp(-par->y) / (par->ss * xP);   // 1 - z_*
    double z = 1.0 - zbar;
    if (z <= 0.0 || zbar <= 0.0) return 0.0;   // this x_P is not possible for this (pD0, y, z_h)

    double qp = pp / z;   // q+_*, Eq. 19
    if (qp > par->qpmax) return 0.0;

    // effective flux: b is not used
    double flux = is_effective_flux(par) ? effective_photon_flux(qp, p) : photon_flux(b, qp, p);

    double xP_dip = min(xP, 0.01);
    double I0, I1;
    exclusive_radial_integrals(pc, m, xP_dip, par->dipole, I0, I1);

    double factor = z*z + zbar*zbar;
    double g = factor * I1*I1 + I0*I0;

    // Jacobian (q+)^2 (1-z) / p+, divided by x_P
    double jac_xpo = qp*qp * zbar / (pp * xP);

    double kernel = z * m*m * jac_xpo * flux * g;

    return kernel * fragmentation(z_h, par);
}

// Exclusive, 5D version: VEGAS integral over {r1, r2, u_qp, b, z_h}. Used at low pt.
double integrand_exclusive_mc(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double r1   = vec[0];
    double r2   = vec[1];
    double u_qp = vec[2];
    double b    = vec[3];
    double z_h   = vec[4];

    double pc = par->pD0 / z_h;   // charm quark's transverse momentum
    double m  = par->m;

    double p2 = pc * pc;
    double mt = sqrt(p2 + par->m2);
    double pp = (mt / M_SQRT2) * exp(par->y);

    // change u_qp (0 to 1) into qp (pp to qpmax)
    double jac_qp = par->qpmax - pp;
    if (jac_qp <= 0.0) return 0.0;
    double qp = pp + u_qp * jac_qp;

    // effective flux: b is not used
    double flux = is_effective_flux(par) ? effective_photon_flux(qp, p) : photon_flux(b, qp, p);

    // light-cone momentum fractions
    double z = pp / qp;
    double zbar = 1.0 - z;

    // W^2, quark-pair mass squared and pomeron momentum fraction
    double w2   = M_SQRT2 * qp * par->ss;
    double Mqq2 = (p2 + par->m2) / (z * zbar);
    double xP   = Mqq2 / w2;
    if (xP <= 0.0 || xP > 0.1) return 0.0;
    double xP_dip = min(xP, 0.01);

    // dipole amplitude for the two dipole sizes
    double D1 = par->dipole->N(r1, xP_dip);
    double D2 = par->dipole->N(r2, xP_dip);

    double f1_r1 = r1 * Jn(1, pc*r1) * Kn(1, m*r1) * D1;
    double f1_r2 = r2 * Jn(1, pc*r2) * Kn(1, m*r2) * D2;
    double f0_r1 = r1 * Jn(0, pc*r1) * Kn(0, m*r1) * D1;
    double f0_r2 = r2 * Jn(0, pc*r2) * Kn(0, m*r2) * D2;

    double factor = z*z + zbar*zbar;
    double g = factor * f1_r1 * f1_r2 + f0_r1 * f0_r2;

    double kernel = z * m*m * jac_qp * flux * g;

    return kernel * fragmentation(z_h, par);
}

// Diffractive: VEGAS integral over {r1, r2, u_k, u_qp, b, u_xpo, z_h}.
double integrand_diffractive(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double k_lo_frag = 0.01;
    double x_po_min_log = 1e-6;

    double r1   = vec[0];
    double r2   = vec[1];
    double k    = k_lo_frag * exp(vec[2] * par->log_k_frag);
    double u_qp = vec[3];
    double b    = vec[4];
    double log_xpo_range = log(0.1 / x_po_min_log);
    double x_po = x_po_min_log * exp(vec[5] * log_xpo_range);
    double z_h   = vec[6];

    double pc = par->pD0 / z_h;
    if (k > pc) return 0.0;   // the gluon cannot have more momentum than the charm quark

    double m2 = par->m2;
    double p2 = pc * pc;
    double mt = sqrt(p2 + m2);
    double pp = (mt / M_SQRT2) * exp(par->y);

    // change u_qp (0 to 1) into qp (pp to qpmax)
    double jac_qp = par->qpmax - pp;
    if (jac_qp <= 0.0) return 0.0;
    double qp = pp + u_qp * jac_qp;

    // effective flux: b is not used
    double flux = is_effective_flux(par) ? effective_photon_flux(qp, p) : photon_flux(b, qp, p);

    // light-cone momentum fractions
    double z = pp / qp;
    double zbar = 1.0 - z;

    // W^2 and quark-pair mass squared
    double w2   = M_SQRT2 * qp * par->ss;
    double Mqq2 = (p2 + m2) / (z * zbar);
    double x    = Mqq2 / (x_po * w2);
    if (x <= 0.0 || x >= 1.0) return 0.0;

    double xfac    = x / (1.0 - x);
    double omega2  = xfac * k * k;
    double sqomega = sqrt(xfac) * k;

    // hard factor H^2
    double denom4      = (p2 + m2) * (p2 + m2);
    double p4m4        = p2*p2 + m2*m2;
    double two_m2p2m2  = 2.0 * m2 * p2;
    double H2 = z / (denom4 * denom4) * ((z*z + zbar*zbar) * p4m4 + two_m2p2m2);

    // gluon diffractive TMD squared
    double gDTMD2 = omega2 * omega2 * r1*r2
                    * Jn(2, k*r1) * Kn(2, sqomega*r1)
                    * Jn(2, k*r2) * Kn(2, sqomega*r2);

    double fun = H2 * (1.0 / (1.0 - x)) * gDTMD2;

    // dipole amplitude for the two dipole sizes
    double x_po_dip = min(x_po, 0.01);
    double D1 = par->dipole->N(r1, x_po_dip);
    double D2 = par->dipole->N(r2, x_po_dip);

    double kernel = 2.0*M_PI * k * jac_qp * flux * D1*D2 * fun
                    * (k * par->log_k_frag) * (x_po * log_xpo_range);

    return kernel * fragmentation(z_h, par);
}

// Diffractive at fixed x_po (par->fixed_xpo):
// VEGAS integral over {r1, r2, u_k, u_qp, b, z_h}.
double integrand_diffractive_xpom(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double k_lo_frag = 0.01;

    double r1   = vec[0];
    double r2   = vec[1];
    double k    = k_lo_frag * exp(vec[2] * par->log_k_frag);
    double u_qp = vec[3];
    double b    = vec[4];
    double z_h   = vec[5];
    double x_po = par->fixed_xpo;

    double pc = par->pD0 / z_h;
    if (k > pc) return 0.0;   

    double m2 = par->m2;
    double p2 = pc * pc;
    double mt = sqrt(p2 + m2);
    double pp = (mt / M_SQRT2) * exp(par->y);

    // change u_qp (0 to 1) into qp (pp to qpmax)
    double jac_qp = par->qpmax - pp;
    if (jac_qp <= 0.0) return 0.0;
    double qp = pp + u_qp * jac_qp;

    // effective flux: b is not used
    double flux = is_effective_flux(par) ? effective_photon_flux(qp, p) : photon_flux(b, qp, p);

    // light-cone momentum fractions
    double z = pp / qp;
    double zbar = 1.0 - z;

    // W^2 and quark-pair mass squared
    double w2   = M_SQRT2 * qp * par->ss;
    double Mqq2 = (p2 + m2) / (z * zbar);
    double x    = Mqq2 / (x_po * w2);
    if (x <= 0.0 || x >= 1.0) return 0.0;

    double xfac    = x / (1.0 - x);
    double omega2  = xfac * k * k;
    double sqomega = sqrt(xfac) * k;

    // hard factor H^2
    double denom4      = (p2 + m2) * (p2 + m2);
    double p4m4        = p2*p2 + m2*m2;
    double two_m2p2m2  = 2.0 * m2 * p2;
    double H2 = z / (denom4 * denom4) * ((z*z + zbar*zbar) * p4m4 + two_m2p2m2);

    // gluon diffractive TMD squared
    double gDTMD2 = omega2 * omega2 * r1*r2
                    * Jn(2, k*r1) * Kn(2, sqomega*r1)
                    * Jn(2, k*r2) * Kn(2, sqomega*r2);

    double fun = H2 * (1.0 / (1.0 - x)) * gDTMD2;

    // dipole amplitude for the two dipole sizes
    double x_po_dip = min(x_po, 0.01);
    double D1 = par->dipole->N(r1, x_po_dip);
    double D2 = par->dipole->N(r2, x_po_dip);

    double kernel = 2.0*M_PI * k * jac_qp * flux * D1*D2 * fun
                    * (k * par->log_k_frag);

    return kernel * fragmentation(z_h, par);
}

// Diffractive at fixed q+ (q+ = 2p+, no flux, no fragmentation):
// VEGAS integral over {r1, r2, u_k, u_xpo}.

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

    double z = 0.5;   // q+ = 2p+, so p+/q+ = 1/2
    double zbar = 1.0 - z;

    double pp = (par->mt / M_SQRT2) * exp(par->y);
    double qp = 2.0 * pp;
    double w2   = M_SQRT2 * qp * par->ss;
    double Mqq2 = (par->p2 + par->m2) / (z * zbar);
    double x    = Mqq2 / (x_po * w2);
    if (x <= 0.0 || x >= 1.0) return 0.0;

    double xfac    = x / (1.0 - x);
    double omega2  = xfac * k * k;
    double sqomega = sqrt(xfac) * k;

    // hard factor H^2
    double H2 = z * par->inv_denom8 * ((z*z + zbar*zbar) * par->p4m4 + par->two_m2p2m2);

    // gluon diffractive TMD squared
    double gDTMD2 = omega2 * omega2 * r1*r2
                    * Jn(2, k*r1) * Kn(2, sqomega*r1)
                    * Jn(2, k*r2) * Kn(2, sqomega*r2);

    double fun = H2 * (1.0 / (1.0 - x)) * gDTMD2;

    // dipole amplitude for the two dipole sizes
    double x_po_dip = min(x_po, 0.01);
    double D1 = par->dipole->N(r1, x_po_dip);
    double D2 = par->dipole->N(r2, x_po_dip);

    double kernel = 2.0*M_PI * k * D1*D2 * fun
                    * (k * par->log_k) * (x_po * log_xpo_range);

    return kernel;
}
