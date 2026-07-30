#ifndef def_hpp
#define def_hpp
#include "amplitudelib.hpp"
#include "interpolation.hpp"
#include <memory>
#include <string>

using namespace std;

enum class FragmentationType { BCFY, KniehlKramer, LHAPDF };

struct parameters
{
    AmplitudeLib *dipole;
    string datafile;
    double m;    // charm quark mass
    double y;
    double ss;
    // Call counts for the Monte Carlo integration (see main.cpp).
    // exclusiveCrossSection (int_exclusive.cpp) picks between two exclusive
    // integrands depending on pD0 vs excl_pt_threshold: below it, the cheap
    // 5D integrand (calls_excl calls -- can afford to be large); at/above
    // it, the factorized 3D one, where r1/r2's oscillatory Bessel integrals
    // are solved analytically per point instead of by VEGAS, so it needs far
    // fewer calls (calls_excl_factorized) but each one costs much more.
    size_t calls_excl;
    size_t calls_excl_factorized;
    double excl_pt_threshold;
    size_t calls_diff;
    double m2;   // just m*m, kept around so we don't recompute it everywhere

    // Photon flux / nuclear parameters
    double alpha, Z, mn, S;
    string channel;
    double bmin, bmax, qpmax;

    // --- Fragmentation (c -> D0 meson), see fragmentation.hpp/lhapdf_grid.hpp ---
    // pD0 is the transverse momentum of the D0 meson we actually observe.
    // The charm quark itself carries a bigger momentum pc = pD0/zh, and we
    // integrate over zh (in [zmin, zmax]) to account for that.
    double pD0;
    double r;              // BCFY non-perturbative parameter
    double N_kk, eps_kk;   // Kniehl & Kramer parameters
    FragmentationType frag_type = FragmentationType::BCFY;
    double zmin, zmax;
    // Pointer to the LHAPDF interpolator (only used when frag_type is LHAPDF).
    // It's just a plain pointer, not a smart pointer, because this whole
    // struct gets copied around a lot and we don't want to own the memory here.
    Interpolator* D_frag_interp = nullptr;

    // Diffractive case only: the emitted gluon's momentum k can't be bigger
    // than the charm quark's own momentum pc = pD0/zh. But pc changes as zh
    // is sampled during the integration, so we can't just set a fixed upper
    // limit for k in the integration box. Instead we make the box for k big
    // enough to cover the largest possible pc (which happens at the smallest
    // zh), and then check "if k > pc, skip this point" inside the integrand.
    double k_upper_frag;
    double log_k_frag;   // log(k_upper_frag / k_lo_frag), computed once before each integration

    // Used only when we fix x_po instead of integrating over it (see
    // integrand_diffractive_xpom / diffractiveCrossSection_xpom). This gives
    // us the cross section at one specific value of x_po.
    double fixed_xpo;

    // Used only in the "fixed W, no photon flux" mode (see main_fixedW.cpp).
    // Here q+ is fixed to 3*p+ instead of being integrated over the photon
    // flux, and there is no fragmentation step -- p is directly the charm
    // quark's transverse momentum K. Since these don't change during the
    // integration, we compute them once outside and just reuse them.
    double p, p2, mt;                             // K, K^2, sqrt(K^2 + m^2)
    double kmax, k_upper, log_k;                  // range for the gluon momentum k
    double denom4, p4m4, two_m2p2m2, inv_denom8;  // pieces used to build H^2 (see integrand_diffractive_fixedW)
};

// D0-level (fragmented) exclusive integrand: dsigma/(d2pD0 dy), evaluated
// at the charm quark momentum pc=par->pD0/zh, with an extra fragmentation
// dimension zh in [par->zmin, par->zmax] and weight D(zh)/zh^2 folded in.
// exclusiveCrossSection (see int_exclusive.cpp) picks between two versions
// depending on pD0:
//  - integrand_exclusive_mc: original 5D VEGAS over {r1, r2, u_qp, b, zh}.
//    Cheap per call, used at low/mid pt.
//  - integrand_exclusive_factorized: r1 and r2 factorize at fixed (u_qp, zh)
//    and are solved with deterministic quadrature inside the integrand
//    instead of being handed to Monte Carlo, leaving a 3D VEGAS over
//    {u_qp, b, zh}. Much better convergence at high pt, where pc=pD0/zh
//    grows large and Jn(pc*r) oscillates too fast for VEGAS to resolve --
//    but each call is more expensive, so it's only worth it there.
double integrand_exclusive_mc(double* vec, size_t dim, void* p);
double integrand_exclusive_factorized(double* vec, size_t dim, void* p);
double exclusiveCrossSection( void* p);

// Diffractive integrand at the D0 level. Also sums over x_po (in [1e-6, 0.1])
// and the gluon momentum k, on top of the fragmentation variable zh.
// 7 integration variables: {r1, r2, u_k, u_qp, b, u_xpo, zh}.
double integrand_diffractive(double* vec, size_t dim, void* p);
double diffractiveCrossSection( void* p);

// Same as integrand_diffractive/diffractiveCrossSection, but x_po is fixed
// (par->fixed_xpo) instead of being one of the integration variables. This
// gives the cross section at one specific x_po, integrated over the
// remaining 6 variables {r1, r2, u_k, u_qp, b, zh}.
double integrand_diffractive_xpom(double* vec, size_t dim, void* p);
double diffractiveCrossSection_xpom(void* p);

double flux_density(double qp, double b, void* p);
double photon_flux(double b, double qp, void* p);

// Bessel function helpers, defined in integrand.cpp. Declared here too so
// int_exclusive.cpp can use them for the fixed-W calculation below.
double Jn(int nu, double x);
double Kn(int nu, double x);

// "Fixed W, no photon flux" mode: q+ is set to 3*p+ by hand instead of being
// integrated over the photon flux, and pc is just par->p directly (no
// fragmentation step). See main_fixedW.cpp for how these are used.
//   exclusive:   computed directly with two 1D integrals (see int_exclusive.cpp),
//                still missing the overall prefactor alpha_em*Nc*e_c^2*sigma0/(2pi^2),
//                same as exclusiveCrossSection above.
//   diffractive: integrated over {r1, r2, u_k, u_xpo}. This one already has
//                the full prefactor alpha_s(mT)*alpha_em*e_c^2*(Nc^2-1)*sigma0/(8pi^4)
//                multiplied in, so nothing else needs to be applied later
//                (except converting GeV^-2 to mb).
double exclusiveCrossSection_fixedW(void* p);
double integrand_diffractive_fixedW(double* vec, size_t dim, void* p);
double diffractiveCrossSection_fixedW(void* p);

#endif
