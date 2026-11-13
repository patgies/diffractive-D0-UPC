#include "def.hpp"
#include "amplitudelib.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>
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


// dN/dq+ = (dN/domega)/sqrt(2)
static double flux_per_qp(double b, double qp, void* p)
{
    parameters* par = (parameters*)p;
    double dN_domega = is_effective_flux(par) ? effective_photon_flux(qp, p) : photon_flux(b, qp, p);
    return dN_domega / M_SQRT2;
}

// D(z_h)/z_h^2
static double fragmentation(double z_h, parameters* par)
{
    double D_frag = par->D_frag_interp->Evaluate(z_h);
    return D_frag / (z_h * z_h);
}

// Charm pc = pD0/z_h and photon q+ in [p+, qpmax] from u_qp in [0, 1]
struct Kinematics {
    double pc, p2;      // pc, pc^2
    double pp, qp;      // p+ of the quark, q+ of the photon
    double jac_qp;
    double z, zbar;     // p+/q+ and 1 - z
    double w2;          // W^2
};

// False if p+ > qpmax
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

// Exclusive r1, r2 integrals as two 1D integrals I0, I1
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

// Exclusive x_P = M_qq^2/W^2; false above XPO_MAX
static bool exclusive_xpo(const Kinematics& kin, const parameters* par, double& xP_dip)
{
    double Mqq2 = (kin.p2 + par->m2) / (kin.z * kin.zbar);
    double xP   = Mqq2 / kin.w2;
    if (xP <= 0.0 || xP > XPO_MAX) return false;
    xP_dip = min(xP, XPO_FROZEN);
    return true;
}

// Exclusive 3D: {u_qp, b, z_h}
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

// Exclusive at fixed x_P: {b, z_h}, returns dsigma/dx_P
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
    if (z <= 0.0 || zbar <= 0.0) return 0.0;   // x_P not reachable

    double qp = pp / z;   // q+ at this x_P
    if (qp > par->qpmax) return 0.0;

    double flux = flux_per_qp(b, qp, p);

    double xP_dip = min(xP, XPO_FROZEN);
    double I0, I1;
    exclusive_radial_integrals(pc, m, xP_dip, par->dipole, I0, I1);

    double factor = z*z + zbar*zbar;
    double g = factor * I1*I1 + I0*I0;

    // Jacobian (q+)^2 (1-z)/p+, over x_P
    double jac_xpo = qp*qp * zbar / (pp * xP);

    double kernel = z * m*m * jac_xpo * flux * g;

    return kernel * fragmentation(z_h, par);
}

// Exclusive 5D: {r1, r2, u_qp, b, z_h}, low pt
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

// Gluon DTMD, arXiv:2606.04169 Eqs. 3.38, 3.43-3.45: 2 pi k H G^2/(1-x) per d ln(x_po),
// G = Omega^2 int dr r J_2(k r) K_2(Omega r) N_g(r), N_g = 2N - N^2, Omega^2 = x k^2/(1-x).
// No k, x_po Jacobians. Returns false outside 0 < x < 1.
static bool gluon_dtmd_kinematics(double k, double x_po, double z, double p2, double w2, const parameters* par,
                                  double& omega, double& pref)
{
    double m2 = par->m2;
    double zbar = 1.0 - z;

    double Mqq2 = (p2 + m2) / (z * zbar);
    double x    = Mqq2 / (x_po * w2);
    if (x <= 0.0 || x >= 1.0) return false;

    double xfac   = x / (1.0 - x);
    double omega2 = xfac * k * k;
    omega = sqrt(xfac) * k;

    // hard factor
    double denom4      = (p2 + m2) * (p2 + m2);
    double p4m4        = p2*p2 + m2*m2;
    double two_m2p2m2  = 2.0 * m2 * p2;
    double H2 = z / (denom4 * denom4) * ((z*z + zbar*zbar) * p4m4 + two_m2p2m2);

    pref = 2.0*M_PI * k * H2 / (1.0 - x) * omega2 * omega2;
    return true;
}

// Gluon DTMD with r1, r2 as VEGAS variables
static double gluon_dtmd_kernel(double r1, double r2, double k, double x_po, double z, double p2,
                                double w2, const parameters* par)
{
    double omega, pref;
    if (!gluon_dtmd_kinematics(k, x_po, z, p2, w2, par, omega, pref)) return 0.0;

    // gluon dipole amplitude
    double x_po_dip = min(x_po, XPO_FROZEN);
    double D1 = par->dipole->N_A(r1, x_po_dip);
    double D2 = par->dipole->N_A(r2, x_po_dip);

    return pref * r1 * Jn(2, k*r1) * Kn(2, omega*r1) * D1
                * r2 * Jn(2, k*r2) * Kn(2, omega*r2) * D2;
}

// Quark DTMD, arXiv:2606.04169 Eqs. 4.21-4.54 at Q^2 = 0 (notes/quark_dtmd.pdf); hard charm + gluon,
// semi-hard charm k: 2 pi k z z(1+zbar^2)/p^4 x/(1-x) (k^2+m^2)/(x k^2+m^2) [Q_1^2 + (1-x) m^2/(x k^2+m^2) Q_0^2],
// Q_n = omega^2 int dr r J_n(k r) K_n(omega r) N(r), omega^2 = (m^2 + x k^2)/(1-x), M^2 = p^2/(z zbar) + m^2/z.
// Per d ln(x_po). Returns false outside 0 < x < 1.
static bool quark_dtmd_kinematics(double k, double x_po, double z, double p2, double w2, const parameters* par,
                                  double& omega, double& c1, double& c0)
{
    double m2 = par->m2;
    double zbar = 1.0 - z;
    double k2 = k * k;

    double M2 = p2 / (z * zbar) + m2 / z;
    double x  = M2 / (x_po * w2);
    if (x <= 0.0 || x >= 1.0) return false;

    double omega2 = (m2 + x * k2) / (1.0 - x);
    omega = sqrt(omega2);

    // hard factor
    double H = z * z * (1.0 + zbar*zbar) / (p2 * p2);

    double xk2m2 = x * k2 + m2;
    c1 = 2.0*M_PI * k * H * x / (1.0 - x) * (k2 + m2) / xk2m2 * omega2 * omega2;
    c0 = c1 * (1.0 - x) * m2 / xk2m2;
    return true;
}

// Quark DTMD with r1, r2 as VEGAS variables
static double quark_dtmd_kernel(double r1, double r2, double k, double x_po, double z, double p2,
                                double w2, const parameters* par)
{
    double omega, c1, c0;
    if (!quark_dtmd_kinematics(k, x_po, z, p2, w2, par, omega, c1, c0)) return 0.0;

    double x_po_dip = min(x_po, XPO_FROZEN);
    double D1 = par->dipole->N(r1, x_po_dip);
    double D2 = par->dipole->N(r2, x_po_dip);
    double radial = r1*r2 * D1*D2;
    double Q1sq = radial * Jn(1, k*r1) * Kn(1, omega*r1) * Jn(1, k*r2) * Kn(1, omega*r2);
    double Q0sq = radial * Jn(0, k*r1) * Kn(0, omega*r1) * Jn(0, k*r2) * Kn(0, omega*r2);

    return c1 * Q1sq + c0 * Q0sq;
}

// k below the cut
static bool below_k_cut(double k, double scale, const parameters* par)
{
    return par->k_cut_factor <= 0.0 || k <= par->k_cut_factor * scale;
}

// Gluon and/or quark DTMD. Gluon cut: k_cut_factor * p or m_T; quark: k_cut_factor * p, p >= quark_p_min.
static double diffractive_kernel(double r1, double r2, double k, double x_po, double z, double p2,
                                 double w2, const parameters* par)
{
    double p = sqrt(p2);
    double gluon_scale = (par->k_cut_scale == "mt") ? sqrt(p2 + par->m2) : p;

    double result = 0.0;
    if (par->dtmd != "quark" && below_k_cut(k, gluon_scale, par))
        result += gluon_dtmd_kernel(r1, r2, k, x_po, z, p2, w2, par);
    if (par->dtmd != "gluon" && p >= par->quark_p_min && below_k_cut(k, p, par))
        result += quark_dtmd_kernel(r1, r2, k, x_po, z, p2, w2, par);
    return result;
}

// Diffractive at given x_po, with flux, fragmentation and Jacobians
static double diffractive_at_xpo(double r1, double r2, double u_k, double u_qp, double b, double z_h,
                                 double x_po, void* p)
{
    parameters* par = (parameters*)p;

    double k = K_MIN * exp(u_k * par->log_k_frag);

    Kinematics kin;
    if (!photon_kinematics(z_h, u_qp, par, kin)) return 0.0;

    double flux = flux_per_qp(b, kin.qp, p);

    double kernel = diffractive_kernel(r1, r2, k, x_po, kin.z, kin.p2, kin.w2, par)
                    * kin.jac_qp * flux * (k * par->log_k_frag);

    return kernel * fragmentation(z_h, par);
}

// Diffractive: {r1, r2, u_k, u_qp, b, u_xpo, z_h}
double integrand_diffractive(double* vec, size_t /*dim*/, void* p)
{
    double log_xpo_range = log(XPO_MAX / XPO_MIN);
    double x_po = XPO_MIN * exp(vec[5] * log_xpo_range);

    return diffractive_at_xpo(vec[0], vec[1], vec[2], vec[3], vec[4], vec[6], x_po, p) * log_xpo_range;
}

// Diffractive at fixed x_po: {r1, r2, u_k, u_qp, b, z_h}
double integrand_diffractive_xpom(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;
    return diffractive_at_xpo(vec[0], vec[1], vec[2], vec[3], vec[4], vec[5], par->fixed_xpo, p);
}

// Fixed q+ = 2p+, no Jacobians
static double diffractive_fixed_qp_kernel(double r1, double r2, double k, double x_po, const parameters* par)
{
    double z = 0.5;   // q+ = 2p+
    double pp = (par->mt / M_SQRT2) * exp(par->y);
    double qp = 2.0 * pp;
    double w2 = M_SQRT2 * qp * par->ss;

    return diffractive_kernel(r1, r2, k, x_po, z, par->p2, w2, par);
}

// Fixed q+: {r1, r2, u_k, u_xpo}
double integrand_diffractive_fixed_qp(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double log_xpo_range = log(XPO_MAX / XPO_MIN);

    double k    = K_MIN * exp(vec[2] * par->log_k);
    double x_po = XPO_MIN * exp(vec[3] * log_xpo_range);

    return diffractive_fixed_qp_kernel(vec[0], vec[1], k, x_po, par)
           * (k * par->log_k) * log_xpo_range;
}

// Fixed q+ and x_po: {r1, r2, u_k}
double integrand_diffractive_fixed_qp_xpom(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double k = K_MIN * exp(vec[2] * par->log_k);

    return diffractive_fixed_qp_kernel(vec[0], vec[1], k, par->fixed_xpo, par) * (k * par->log_k);
}

// int dr r J_nu(k r) K_nu(omega r) N(r) (N_A if adjoint), Gauss-Legendre panels on a geometric grid
// plus the half periods pi/k. rmax = min(99, 50/omega).
double dtmd_radial_integral(int nu, double k, double omega, double x_dip, bool adjoint, AmplitudeLib* dipole)
{
    static const gsl_integration_glfixed_table* gl = gsl_integration_glfixed_table_alloc(10);
    static const int refine = getenv("QUAD_REFINE") ? atoi(getenv("QUAD_REFINE")) : 1;
    const double rmin = 1e-4;
    const int n_geometric = 40 * refine;
    double rmax = min(99.0, 50.0 / omega);

    vector<double> edges = {0.0};
    for (int i = 0; i <= n_geometric; ++i)
        edges.push_back(rmin * pow(rmax / rmin, (double)i / n_geometric));
    for (double r = M_PI / (k * refine); r < rmax; r += M_PI / (k * refine))
        edges.push_back(r);
    sort(edges.begin(), edges.end());

    double sum = 0.0;
    for (size_t i = 0; i + 1 < edges.size(); ++i) {
        double a = edges[i], b = edges[i + 1];
        if (b - a < 1e-12) continue;
        for (size_t j = 0; j < gl->n; ++j) {
            double r, w;
            gsl_integration_glfixed_point(a, b, j, &r, &w, gl);
            double N = adjoint ? dipole->N_A(r, x_dip) : dipole->N(r, x_dip);
            sum += w * r * Jn(nu, k*r) * Kn(nu, omega*r) * N;
        }
    }
    return sum;
}

// int_a^b f, 8-point Gauss-Legendre panels
template <class F>
static double gauss_legendre(double a, double b, int n_panels, F f)
{
    static const gsl_integration_glfixed_table* gl = gsl_integration_glfixed_table_alloc(8);
    double sum = 0.0, h = (b - a) / n_panels;
    for (int i = 0; i < n_panels; ++i)
        for (size_t j = 0; j < gl->n; ++j) {
            double t, w;
            gsl_integration_glfixed_point(a + i*h, a + (i+1)*h, j, &t, &w, gl);
            sum += w * f(t);
        }
    return sum;
}

// Fixed q+ by quadrature (no Monte Carlo); same as diffractiveCrossSection_fixed_qp
double diffractiveCrossSection_fixed_qp_quad(void* p)
{
    parameters* par = (parameters*)p;
    // QUAD_REFINE: more panels, for convergence checks
    const int refine = getenv("QUAD_REFINE") ? atoi(getenv("QUAD_REFINE")) : 1;
    const int n_panels_k = 16 * refine, n_panels_xpo = 12 * refine;

    double z = 0.5;   // q+ = 2p+
    double zbar = 1.0 - z;
    double pp = (par->mt / M_SQRT2) * exp(par->y);
    double w2 = M_SQRT2 * 2.0 * pp * par->ss;
    double pc = par->p, p2 = par->p2;

    // int d ln(x_po) d ln(k) k f, k up to the cut, x_po from M^2/W^2
    auto integrate = [&](double scale, double M2, auto f) {
        double k_hi = (par->k_cut_factor > 0.0) ? min(par->k_upper, par->k_cut_factor * scale) : par->k_upper;
        double xpo_lo = max(XPO_MIN, M2 / w2);
        if (k_hi <= K_MIN || xpo_lo >= XPO_MAX) return 0.0;
        return gauss_legendre(log(xpo_lo), log(XPO_MAX), n_panels_xpo, [&](double ln_xpo) {
            double x_po = exp(ln_xpo);
            return gauss_legendre(log(K_MIN), log(k_hi), n_panels_k, [&](double ln_k) {
                double k = exp(ln_k);
                return k * f(k, x_po);
            });
        });
    };

    double result = 0.0;
    if (par->dtmd != "quark") {
        double gluon_scale = (par->k_cut_scale == "mt") ? par->mt : pc;
        result += integrate(gluon_scale, (p2 + par->m2) / (z * zbar), [&](double k, double x_po) {
            double omega, pref;
            if (!gluon_dtmd_kinematics(k, x_po, z, p2, w2, par, omega, pref)) return 0.0;
            double I = dtmd_radial_integral(2, k, omega, min(x_po, XPO_FROZEN), true, par->dipole);
            return pref * I*I;
        });
    }
    if (par->dtmd != "gluon" && pc >= par->quark_p_min) {
        result += integrate(pc, p2 / (z * zbar) + par->m2 / z, [&](double k, double x_po) {
            double omega, c1, c0;
            if (!quark_dtmd_kinematics(k, x_po, z, p2, w2, par, omega, c1, c0)) return 0.0;
            double x_dip = min(x_po, XPO_FROZEN);
            double I1 = dtmd_radial_integral(1, k, omega, x_dip, false, par->dipole);
            double I0 = dtmd_radial_integral(0, k, omega, x_dip, false, par->dipole);
            return c1 * I1*I1 + c0 * I0*I0;
        });
    }
    return result;
}
