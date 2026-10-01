#ifndef def_hpp
#define def_hpp
#include "amplitudelib.hpp"
#include "interpolation.hpp"
#include <memory>
#include <string>

using namespace std;

enum class FragmentationType { BCFY, KniehlKramer, HymnD };

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
    // "AA" (default): Pb-Pb, Gamma_AA(b) tabulated survival probability +
    // EMD channel factors (An0n/Xn0n), as everywhere above. "pA": proton
    // target (Sec. 6 of arXiv:2606.05469) -- the SAME lead-ion flux, but
    // Gamma_pA(b)=exp(-sigma_NN*T_A(b)) in its place (see init_pA_flux() /
    // GammaPA() in photon_flux.cpp) and no EMD factor at all. Only affects
    // photon_flux(); not meaningful together with flux_model=="EFF" (that
    // convolution is over a nuclear target's own extent, which a proton
    // target doesn't have).
    string target = "AA";
    // "PL" (default): point-like bare flux, K1(eta)^2 (flux_density() in
    // photon_flux.cpp). "WS": Woods-Saxon bare flux (flux_density_WS()), the
    // Vidovic-Greiner-Best-Soff / Krauss-Greiner-Soff form with the nuclear
    // charge form factor folded in -- see photon_flux.hpp for the derivation.
    string flux_model = "PL";

    // --- Fragmentation (c -> D0 meson), see fragmentation.hpp/hymnd_grid.hpp ---
    // pD0 is the transverse momentum of the D0 meson we actually observe.
    // The charm quark itself carries a bigger momentum pc = pD0/zh, and we
    // integrate over zh (in [zmin, zmax]) to account for that.
    double pD0;
    double r;              // BCFY non-perturbative parameter
    double N_kk, eps_kk;   // Kniehl & Kramer parameters
    FragmentationType frag_type = FragmentationType::BCFY;
    double zmin, zmax;
    // Pointer to the HymnD interpolator (only used when frag_type is HymnD).
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

    // Used only in the "fixed q+, no photon flux" mode (see main_fixed_qp.cpp).
    // Here q+ is fixed (2*p+ in D0_fixed_qp, 3*p+ in D0_fixed_qp_xpom) instead of the photon
    // flux, and there is no fragmentation step -- p is directly the charm
    // quark's transverse momentum K. Since these don't change during the
    // integration, we compute them once outside and just reuse them.
    double p, p2, mt;                             // K, K^2, sqrt(K^2 + m^2)
    double kmax, k_upper, log_k;                  // range for the gluon momentum k
    double denom4, p4m4, two_m2p2m2, inv_denom8;  // pieces used to build H^2 (see integrand_diffractive_fixed_qp)
};

// Lower edge of the b integration. Only the point-like PL(AnAn) channel
// (and the pA "PL" comparison) has a sharp cut at bmin, since there
// Gamma = 1; every channel with a hadronic survival factor Gamma_AA(b)
// starts at b = 0 -- Gamma already removes the overlapping configurations.
// "EFF"/"STARLIGHT": b is an unused placeholder in [0,1].
inline bool is_effective_flux(const parameters* par)
{
    // b already integrated out: our own effective flux ("EFF") or
    // Paakkinen's tabulated one ("STARLIGHT", see init_starlight_flux())
    return par->flux_model == "EFF" || par->flux_model == "STARLIGHT";
}

inline double b_integration_min(const parameters* par)
{
    if (is_effective_flux(par)) return 0.0;
    if (par->target == "pA" || par->channel == "PL(AnAn)") return par->bmin;
    return 0.0;
}

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

// Same as exclusiveCrossSection, but at one fixed x_P (par->fixed_xpo)
// instead of integrating over the photon energy q+ -- x_P isn't an
// independent variable for the exclusive channel (unlike the diffractive
// channel's x_po), so fixing it picks out one specific q+ via a delta
// function (see integrand_exclusive_xpom), leaving {b, zh} for VEGAS.
double integrand_exclusive_xpom(double* vec, size_t dim, void* p);
double exclusiveCrossSection_xpom(void* p);

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
double flux_density_WS(double z, double b, void* p);
double photon_flux(double b, double qp, void* p);
// "Effective flux": drop-in replacement for \int db photon_flux(b,qp,par) db
// under the geometric convolution of arXiv:2404.09731 Eq. 4 (target-nucleus
// spatial extent folded in) instead of the single-b treatment above --
// see photon_flux.hpp/init_effective_flux() for how it's built. Takes no b:
// the whole b-dependence is already integrated out.
double effective_photon_flux(double qp, void* p);

// Bessel function helpers, defined in integrand.cpp. Declared here too so
// int_exclusive.cpp can use them for the fixed-q+ calculation below.
double Jn(int nu, double x);
double Kn(int nu, double x);

// "Fixed q+, no photon flux" mode: q+ is set by hand (2*p+, or 3*p+ for _xpom) instead of being
// integrated over the photon flux, and pc is just par->p directly (no
// fragmentation step). See main_fixed_qp.cpp for how these are used.
//   exclusive:   computed directly with two 1D integrals (see int_exclusive.cpp),
//                still missing the overall prefactor alpha_em*Nc*e_c^2*sigma0/(2pi^2),
//                same as exclusiveCrossSection above.
//   diffractive: integrated over {r1, r2, u_k, u_xpo}. This one already has
//                the full prefactor alpha_s(mT)*alpha_em*e_c^2*(Nc^2-1)*sigma0/(8pi^4)
//                multiplied in, so nothing else needs to be applied later
//                (except converting GeV^-2 to mb).
double exclusiveCrossSection_fixed_qp(void* p);
double integrand_diffractive_fixed_qp(double* vec, size_t dim, void* p);
double diffractiveCrossSection_fixed_qp(void* p);

// Same as the fixed_qp pair above, but x_po is fixed at par->fixed_xpo
// instead of being integrated over (dropping u_xpo), giving
// dsigma_fixed_qp/(d2K dx_po) at one specific x_po -- the fixed-q+ analogue
// of diffractiveCrossSection_xpom.
double integrand_diffractive_fixed_qp_xpom(double* vec, size_t dim, void* p);
double diffractiveCrossSection_fixed_qp_xpom(void* p);

#endif
