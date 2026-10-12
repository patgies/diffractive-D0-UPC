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


// Photon flux per unit q+. The flux functions return dN/domega and q+ = sqrt(2) omega, so
// dN/dq+ = (dN/domega)/sqrt(2).
static double flux_per_qp(double b, double qp, void* p)
{
    parameters* par = (parameters*)p;
    double dN_domega = is_effective_flux(par) ? effective_photon_flux(qp, p) : photon_flux(b, qp, p);
    return dN_domega / M_SQRT2;
}

// Fragmentation factor D(z_h)/z_h^2 for c -> D0.
static double fragmentation(double z_h, parameters* par)
{
    double D_frag = par->D_frag_interp->Evaluate(z_h);
    return D_frag / (z_h * z_h);
}

// Charm quark with transverse momentum pc = pD0/z_h and rapidity par->y, and a photon with momentum qp
// between p+ and qpmax (u_qp from 0 to 1). jac_qp is the Jacobian of u_qp.
struct Kinematics {
    double pc, p2;      // quark transverse momentum and its square
    double pp, qp;      // p+ of the quark, q+ of the photon
    double jac_qp;
    double z, zbar;     // p+/q+ and 1 - z
    double w2;          // W^2
};

// Returns false if no photon energy is possible (p+ above qpmax).
static bool photon_kinematics(double z_h, double u_qp, const parameters* par, Kinematics& kin)
{
    kin.pc = par->pD0 / z_h;
    kin.p2 = kin.pc * kin.pc;
    double mt = sqrt(kin.p2 + par->m2);
    kin.pp = (mt / M_SQRT2) * exp(par->y);

    kin.jac_qp = par->qpmax - kin.pp;
    if (kin.jac_qp <= 0.0) return false;
    kin.qp = kin.pp + u_qp * kin.jac_qp;

    kin.z = kin.pp / kin.qp;
    kin.zbar = 1.0 - kin.z;
    kin.w2 = M_SQRT2 * kin.qp * par->ss;
    return true;
}

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

// The r1 and r2 integrals of the exclusive process as two 1D integrals, I0 and I1.
void exclusive_radial_integrals(double pc, double m, double x_dip, AmplitudeLib* dipole, double& I0, double& I1)
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

// x_P of the exclusive process, M_qq^2/W^2. Returns false above XPO_MAX.
static bool exclusive_xpo(const Kinematics& kin, const parameters* par, double& xP_dip)
{
    double Mqq2 = (kin.p2 + par->m2) / (kin.z * kin.zbar);
    double xP   = Mqq2 / kin.w2;
    if (xP <= 0.0 || xP > XPO_MAX) return false;
    xP_dip = min(xP, XPO_FROZEN);
    return true;
}

// Exclusive, 3D version: VEGAS integral over {u_qp, b, z_h}. r1 and r2 are done inside.
double integrand_exclusive_factorized(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double u_qp = vec[0];
    double b    = vec[1];
    double z_h  = vec[2];

    Kinematics kin;
    if (!photon_kinematics(z_h, u_qp, par, kin)) return 0.0;

    double flux = flux_per_qp(b, kin.qp, p);

    double xP_dip;
    if (!exclusive_xpo(kin, par, xP_dip)) return 0.0;

    double I0, I1;
    exclusive_radial_integrals(kin.pc, par->m, xP_dip, par->dipole, I0, I1);

    double factor = kin.z*kin.z + kin.zbar*kin.zbar;
    double g = factor * I1*I1 + I0*I0;

    double kernel = kin.z * par->m*par->m * kin.jac_qp * flux * g;

    return kernel * fragmentation(z_h, par);
}

// Exclusive at fixed x_P (par->fixed_xpo): x_P fixes q+, so VEGAS only
// integrates over {b, z_h}. Returns dsigma/dx_P.
double integrand_exclusive_xpom(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double b   = vec[0];
    double z_h = vec[1];

    double xP = par->fixed_xpo;
    if (xP <= 0.0 || xP > XPO_MAX) return 0.0;

    double pc = par->pD0 / z_h;
    double m  = par->m;
    double p2 = pc * pc;
    double mt = sqrt(p2 + par->m2);
    double pp = (mt / M_SQRT2) * exp(par->y);

    double zbar = mt * exp(-par->y) / (par->ss * xP);   // 1 - z_*
    double z = 1.0 - zbar;
    if (z <= 0.0 || zbar <= 0.0) return 0.0;   // this x_P is not possible for this (pD0, y, z_h)

    double qp = pp / z;   // q+ at this x_P
    if (qp > par->qpmax) return 0.0;

    double flux = flux_per_qp(b, qp, p);

    double xP_dip = min(xP, XPO_FROZEN);
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
    double z_h  = vec[4];

    Kinematics kin;
    if (!photon_kinematics(z_h, u_qp, par, kin)) return 0.0;

    double flux = flux_per_qp(b, kin.qp, p);

    double xP_dip;
    if (!exclusive_xpo(kin, par, xP_dip)) return 0.0;

    double D1 = par->dipole->N(r1, xP_dip);
    double D2 = par->dipole->N(r2, xP_dip);

    double pc = kin.pc, m = par->m;
    double f1_r1 = r1 * Jn(1, pc*r1) * Kn(1, m*r1) * D1;
    double f1_r2 = r2 * Jn(1, pc*r2) * Kn(1, m*r2) * D2;
    double f0_r1 = r1 * Jn(0, pc*r1) * Kn(0, m*r1) * D1;
    double f0_r2 = r2 * Jn(0, pc*r2) * Kn(0, m*r2) * D2;

    double factor = kin.z*kin.z + kin.zbar*kin.zbar;
    double g = factor * f1_r1 * f1_r2 + f0_r1 * f0_r2;

    double kernel = kin.z * m*m * kin.jac_qp * flux * g;

    return kernel * fragmentation(z_h, par);
}

// Diffractive TMD cross section at given photon energy (W^2 = w2), quark momentum (p2 = p^2) and
// z = p+/q+, as in arXiv:2606.04169 (Eqs. 3.38, 3.43-3.45): 2 pi k H(z, p) G(r1) G(r2) / (1 - x),
// where G has the gluon-gluon dipole amplitude N_g = 2N - N^2 (N_A) and x = M_qq^2/(x_po W^2). The
// cross section is per d ln(1/x) = d ln(x_po). The Jacobians of the k and x_po variables are not included.
static double diffractive_kernel(double r1, double r2, double k, double x_po, double z, double p2,
                                 double w2, const parameters* par)
{
    double m2 = par->m2;
    double zbar = 1.0 - z;

    double Mqq2 = (p2 + m2) / (z * zbar);
    double x    = Mqq2 / (x_po * w2);
    if (x <= 0.0 || x >= 1.0) return 0.0;

    double xfac    = x / (1.0 - x);
    double omega2  = xfac * k * k;
    double sqomega = sqrt(xfac) * k;

    // hard factor H
    double denom4      = (p2 + m2) * (p2 + m2);
    double p4m4        = p2*p2 + m2*m2;
    double two_m2p2m2  = 2.0 * m2 * p2;
    double H2 = z / (denom4 * denom4) * ((z*z + zbar*zbar) * p4m4 + two_m2p2m2);

    // gluon diffractive TMD
    double gDTMD2 = omega2 * omega2 * r1*r2
                    * Jn(2, k*r1) * Kn(2, sqomega*r1)
                    * Jn(2, k*r2) * Kn(2, sqomega*r2);

    double fun = H2 * (1.0 / (1.0 - x)) * gDTMD2;

    // gluon dipole amplitude
    double x_po_dip = min(x_po, XPO_FROZEN);
    double D1 = par->dipole->N_A(r1, x_po_dip);
    double D2 = par->dipole->N_A(r2, x_po_dip);

    return 2.0*M_PI * k * D1*D2 * fun;
}

// Diffractive at given x_po, for the VEGAS variables {r1, r2, u_k, u_qp, b, z_h}: the kernel with the
// photon flux, the fragmentation function and the Jacobians of u_qp and u_k.
static double diffractive_at_xpo(double r1, double r2, double u_k, double u_qp, double b, double z_h,
                                 double x_po, void* p)
{
    parameters* par = (parameters*)p;

    double k = K_MIN * exp(u_k * par->log_k_frag);

    Kinematics kin;
    if (k > par->pD0 / z_h) return 0.0;   // the gluon cannot have more momentum than the charm quark
    if (!photon_kinematics(z_h, u_qp, par, kin)) return 0.0;

    double flux = flux_per_qp(b, kin.qp, p);

    double kernel = diffractive_kernel(r1, r2, k, x_po, kin.z, kin.p2, kin.w2, par)
                    * kin.jac_qp * flux * (k * par->log_k_frag);

    return kernel * fragmentation(z_h, par);
}

// Diffractive: VEGAS integral over {r1, r2, u_k, u_qp, b, u_xpo, z_h}.
double integrand_diffractive(double* vec, size_t /*dim*/, void* p)
{
    double log_xpo_range = log(XPO_MAX / XPO_MIN);
    double x_po = XPO_MIN * exp(vec[5] * log_xpo_range);

    return diffractive_at_xpo(vec[0], vec[1], vec[2], vec[3], vec[4], vec[6], x_po, p) * log_xpo_range;
}

// Diffractive at fixed x_po (par->fixed_xpo): VEGAS integral over {r1, r2, u_k, u_qp, b, z_h}.
double integrand_diffractive_xpom(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;
    return diffractive_at_xpo(vec[0], vec[1], vec[2], vec[3], vec[4], vec[5], par->fixed_xpo, p);
}

// Diffractive at fixed q+ (q+ = 2p+, no flux, no fragmentation), at given r1, r2, k and x_po,
// without the Jacobians of the k and x_po variables.
static double diffractive_fixed_qp_kernel(double r1, double r2, double k, double x_po, const parameters* par)
{
    if (k > par->p) return 0.0;

    double z = 0.5;   // q+ = 2p+, so p+/q+ = 1/2
    double pp = (par->mt / M_SQRT2) * exp(par->y);
    double qp = 2.0 * pp;
    double w2 = M_SQRT2 * qp * par->ss;

    return diffractive_kernel(r1, r2, k, x_po, z, par->p2, w2, par);
}

// Diffractive at fixed q+: VEGAS integral over {r1, r2, u_k, u_xpo}.
double integrand_diffractive_fixed_qp(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double log_xpo_range = log(XPO_MAX / XPO_MIN);

    double k    = K_MIN * exp(vec[2] * par->log_k);
    double x_po = XPO_MIN * exp(vec[3] * log_xpo_range);

    return diffractive_fixed_qp_kernel(vec[0], vec[1], k, x_po, par)
           * (k * par->log_k) * log_xpo_range;
}

// Diffractive at fixed q+ and fixed x_po (par->fixed_xpo): VEGAS integral over {r1, r2, u_k}.
double integrand_diffractive_fixed_qp_xpom(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double k = K_MIN * exp(vec[2] * par->log_k);

    return diffractive_fixed_qp_kernel(vec[0], vec[1], k, par->fixed_xpo, par) * (k * par->log_k);
}
