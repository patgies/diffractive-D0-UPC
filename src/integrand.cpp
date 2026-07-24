#include "def.hpp"
#include "amplitudelib.hpp"
#include "fragmentation.hpp"
#include <cmath>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_bessel.h>

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


// D(zh)/zh^2 fragmentation (c -> D0)
static double fragmentation(double zh, parameters* par)
{
    double D_frag;
    switch (par->frag_type) {
        case FragmentationType::KniehlKramer:
            D_frag = D_kniehl_kramer(zh, par->N_kk, par->eps_kk);
            break;
        case FragmentationType::LHAPDF:
            D_frag = par->D_frag_interp->Evaluate(zh);
            break;
        case FragmentationType::BCFY:
        default:
            D_frag = Dc_to_D0(zh, par->r);
            break;
    }
    return D_frag / (zh * zh);
}

// D0-level exclusive integrand: dsigma/(d2pD0 dy)

double integrand_exclusive(double* vec, size_t /*dim*/, void* p)
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

    // qp from unit interval, in [pp, qpmax]
    double jac_qp = par->qpmax - pp;
    if (jac_qp <= 0.0) return 0.0;
    double qp = pp + u_qp * jac_qp;

    // Photon flux
    double flux = photon_flux(b, qp, p);

    // Light-cone fractions
    double z1 = pp / qp;
    double z2 = 1.0 - z1;

    // W^2, Mqq^2, xP
    double w2   = M_SQRT2 * qp * par->ss;
    double Mqq2 = (p2 + par->m2) / (z1 * z2);
    double xP   = Mqq2 / w2;
    if (xP <= 0.0 || xP > 0.1) return 0.0;
    double xP_dip = std::min(xP, 0.01);

    // Dipole amplitudes 
    double D1 = par->dipole->N(r1, xP_dip);
    double D2 = par->dipole->N(r2, xP_dip);

    // Radial kernels
    double f1_r1 = r1 * Jn(1, pc*r1) * Kn(1, m*r1) * D1;
    double f1_r2 = r2 * Jn(1, pc*r2) * Kn(1, m*r2) * D2;
    double f0_r1 = r1 * Jn(0, pc*r1) * Kn(0, m*r1) * D1;
    double f0_r2 = r2 * Jn(0, pc*r2) * Kn(0, m*r2) * D2;

    double factor = z1*z1 + z2*z2;
    double g = factor * f1_r1 * f1_r2 + f0_r1 * f0_r2;

    double kernel = z1 * m*m * jac_qp * flux * g;

    return kernel * fragmentation(zh, par);
}

// D0-level (fragmented) diffractive integrand

double integrand_diffractive(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    constexpr double k_lo_frag = 0.01;    
    constexpr double x_po_min_log = 1e-6;

    double r1   = vec[0];
    double r2   = vec[1];
    double k    = k_lo_frag * exp(vec[2] * par->log_k_frag);
    double u_qp = vec[3];
    double b    = vec[4];
    double log_xpo_range = log(0.1 / x_po_min_log);
    double x_po = x_po_min_log * exp(vec[5] * log_xpo_range);
    double zh   = vec[6];

    double pc = par->pD0 / zh;
    if (k > pc) return 0.0;   // gluon can't carry more transverse momentum than the charm quark

    double m2 = par->m2;
    double p2 = pc * pc;
    double mt = sqrt(p2 + m2);
    double pp = (mt / M_SQRT2) * exp(par->y);

    // qp from unit interval
    double jac_qp = par->qpmax - pp;
    if (jac_qp <= 0.0) return 0.0;
    double qp = pp + u_qp * jac_qp;

    // Photon flux
    double flux = photon_flux(b, qp, p);

    // Light-cone fractions
    double z1 = pp / qp;
    double z2 = 1.0 - z1;

    // W^2, Mqq^2, x
    double w2   = M_SQRT2 * qp * par->ss;
    double Mqq2 = (p2 + m2) / (z1 * z2);
    double x    = Mqq2 / (x_po * w2);
    if (x <= 0.0 || x >= 1.0) return 0.0;

    // Transverse scale
    double xfac    = x / (1.0 - x);
    double omega2  = xfac * k * k;
    double sqomega = sqrt(xfac) * k;

    // Hard factor H^2 
    double denom4      = (p2 + m2) * (p2 + m2);
    double p4m4        = p2*p2 + m2*m2;
    double two_m2p2m2  = 2.0 * m2 * p2;
    double H2 = z1 / (denom4 * denom4) * ((z1*z1 + z2*z2) * p4m4 + two_m2p2m2);

    // gDTMD^2 
    double gDTMD2 = omega2 * omega2 * r1*r2
                    * Jn(2, k*r1) * Kn(2, sqomega*r1)
                    * Jn(2, k*r2) * Kn(2, sqomega*r2);

    double fun = H2 * (1.0 / (1.0 - x)) * gDTMD2;

    // Dipole amplitudes
    double x_po_dip = std::min(x_po, 0.01);
    double D1 = par->dipole->N(r1, x_po_dip);
    double D2 = par->dipole->N(r2, x_po_dip);

    double kernel = 2.0*M_PI * k * jac_qp * flux * D1*D2 * fun
                    * (k * par->log_k_frag) * (x_po * log_xpo_range);

    return kernel * fragmentation(zh, par);
}
