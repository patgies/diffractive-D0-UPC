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
    // 3D integral at or above it (see cross_section_exclusive.cpp).
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
    // has pc = pD0/z_h, and we integrate over z_h from z_h_min to z_h_max.
    double pD0;
    FragmentationType frag_type = FragmentationType::BCFY;
    double z_h_min, z_h_max;
    // D(z_h) at the fragmentation scale (the memory belongs to main).
    Interpolator* D_frag_interp = nullptr;

    // Diffractive: the gluon momentum k goes up to k_upper_frag = pD0/z_h_min,
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

// Exclusive, dsigma/(d2pD0 dy): 5D VEGAS integral over {r1, r2, u_qp, b, z_h}
// at low pt, 3D integral over {u_qp, b, z_h} at high pt.
double integrand_exclusive_mc(double* vec, size_t dim, void* p);
double integrand_exclusive_factorized(double* vec, size_t dim, void* p);
double exclusiveCrossSection( void* p);

// Exclusive at fixed x_P (par->fixed_xpo): VEGAS integral over {b, z_h}.
double integrand_exclusive_xpom(double* vec, size_t dim, void* p);
double exclusiveCrossSection_xpom(void* p);

// Diffractive: 7D VEGAS integral over {r1, r2, u_k, u_qp, b, u_xpo, z_h}.
double integrand_diffractive(double* vec, size_t dim, void* p);
double diffractiveCrossSection( void* p);

// Diffractive at fixed x_po: 6D VEGAS integral over {r1, r2, u_k, u_qp, b, z_h}.
double integrand_diffractive_xpom(double* vec, size_t dim, void* p);
double diffractiveCrossSection_xpom(void* p);

// Impact parameters of the flux (as in arXiv:2404.09731):
//   r = distance from the photon to the centre of the nucleus that emits it,
//   b = distance between the centres of the two nuclei (Gamma_AA and the EMD factor depend on it),
//   s = position in the target nucleus (only in the effective flux, b = |r - s|).
// The impact parameter of the dipole amplitude, b_d, is not a variable of the C++ code:
// it is fixed by the Glauber file glauber_mve_<b_d>, and integrated in plotting_scripts/.
double flux_density(double qp, double r, void* p);
double flux_density_WS(double z_gamma, double r, void* p);
double photon_flux(double b, double qp, void* p);
// Effective flux (arXiv:2404.09731 Eq. 4): used in place of the b integral of photon_flux().
double effective_photon_flux(double qp, void* p);

// Bessel functions, defined in integrands.cpp.
double Jn(int nu, double x);
double Kn(int nu, double x);

// Fixed-q+ mode (q+ = 2p+, no flux, no fragmentation). The results do not have
// the prefactors yet (they are applied in plotting_scripts/charm_fixed_qp.py).
double exclusiveCrossSection_fixed_qp(void* p);
double integrand_diffractive_fixed_qp(double* vec, size_t dim, void* p);
double diffractiveCrossSection_fixed_qp(void* p);

#endif
