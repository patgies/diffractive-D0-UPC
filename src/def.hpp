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
    double m;    // charm quark mass
    double y;
    double ss;
    // Number of VEGAS calls. Exclusive: 5D integral below excl_pt_threshold,
    // 3D integral at or above it (see int_exclusive.cpp).
    size_t calls_excl;
    size_t calls_excl_factorized;
    double excl_pt_threshold;
    size_t calls_diff;
    double m2;   // m*m, saved so it is not computed again each time

    // Photon flux and nucleus parameters
    double alpha, Z, mn, S;
    string channel;
    double bmin, bmax, qpmax;
    // "AA" (default): Pb-Pb, uses Gamma_AA(b) and the EMD factor of the channel.
    // "pA": proton target, uses Gamma_pA(b) and no EMD factor.
    string target = "AA";
    // Flux model: "PL" point-like, "WS" Woods-Saxon. "EFF" and "STARLIGHT" are
    // effective fluxes: the sum over b is already done.
    string flux_model = "PL";

    // Fragmentation c -> D0: pD0 is the D0 momentum we measure. The charm quark
    // has pc = pD0/zh, and we integrate over zh from zmin to zmax.
    double pD0;
    FragmentationType frag_type = FragmentationType::BCFY;
    double zmin, zmax;
    // D(z) at the fragmentation scale (the memory belongs to main).
    Interpolator* D_frag_interp = nullptr;

    // Diffractive: the gluon momentum k goes up to k_upper_frag = pD0/zmin,
    // and points with k > pc give zero.
    double k_upper_frag;
    double log_k_frag;   // log(k_upper_frag / k_lo_frag), set before each integration

    // Fixed value of x_po, used by the *_xpom functions.
    double fixed_xpo;

    // Fixed-q+ mode (main_fixed_qp.cpp): no flux, no fragmentation;
    // p is the charm quark's transverse momentum K.
    double p, p2, mt;                             // K, K^2, sqrt(K^2 + m^2)
    double kmax, k_upper, log_k;                  // range for the gluon momentum k
    double denom4, p4m4, two_m2p2m2, inv_denom8;  // parts of H^2 (see integrand_diffractive_fixed_qp)
};

// For "EFF"/"STARLIGHT" the sum over b is already done. If not, the b integral starts
// at bmin for PL(AnAn) and pA "PL", and at 0 for the other channels.
inline bool is_effective_flux(const parameters* par)
{
    return par->flux_model == "EFF" || par->flux_model == "STARLIGHT";
}

inline double b_integration_min(const parameters* par)
{
    if (is_effective_flux(par)) return 0.0;
    if (par->target == "pA" || par->channel == "PL(AnAn)") return par->bmin;
    return 0.0;
}

// Exclusive, dsigma/(d2pD0 dy): 5D VEGAS integral over {r1, r2, u_qp, b, zh}
// at low pt, 3D integral over {u_qp, b, zh} at high pt.
double integrand_exclusive_mc(double* vec, size_t dim, void* p);
double integrand_exclusive_factorized(double* vec, size_t dim, void* p);
double exclusiveCrossSection( void* p);

// Exclusive at fixed x_P (par->fixed_xpo): VEGAS integral over {b, zh}.
double integrand_exclusive_xpom(double* vec, size_t dim, void* p);
double exclusiveCrossSection_xpom(void* p);

// Diffractive: 7D VEGAS integral over {r1, r2, u_k, u_qp, b, u_xpo, zh}.
double integrand_diffractive(double* vec, size_t dim, void* p);
double diffractiveCrossSection( void* p);

// Diffractive at fixed x_po: 6D VEGAS integral over {r1, r2, u_k, u_qp, b, zh}.
double integrand_diffractive_xpom(double* vec, size_t dim, void* p);
double diffractiveCrossSection_xpom(void* p);

double flux_density(double qp, double b, void* p);
double flux_density_WS(double z, double b, void* p);
double photon_flux(double b, double qp, void* p);
// Effective flux (arXiv:2404.09731 Eq. 4): used in place of the b integral of photon_flux().
double effective_photon_flux(double qp, void* p);

// Bessel functions, defined in integrand.cpp.
double Jn(int nu, double x);
double Kn(int nu, double x);

// Fixed-q+ mode (q+ = 2p+, no flux, no fragmentation). The exclusive result
// does not have the prefactor yet. The diffractive result has it.
double exclusiveCrossSection_fixed_qp(void* p);
double integrand_diffractive_fixed_qp(double* vec, size_t dim, void* p);
double diffractiveCrossSection_fixed_qp(void* p);

#endif
