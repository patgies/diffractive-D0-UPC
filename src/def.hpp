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
    double m;
    double y;
    double ss;
    size_t calls;
    double m2;

    // Photon flux / nuclear parameters
    double alpha, Z, mn, S;
    string channel;
    double bmin, bmax, qpmax;

    // --- Fragmentation (c -> D0), see fragmentation.hpp/lhapdf_grid.hpp ---
    // pD0 is the *observed D0 meson's* transverse momentum (the grid input
    // for the _frag integrands below); the charm quark that actually
    // scatters carries pc = pD0/zh, a different, larger momentum integrated
    // out via the fragmentation variable zh in [zmin, zmax].
    double pD0;
    double r;              // BCFY non-perturbative parameter
    double N_kk, eps_kk;   // Kniehl & Kramer parameters
    FragmentationType frag_type = FragmentationType::BCFY;
    double zmin, zmax;
    // Raw, non-owning pointer (same convention as `dipole` above): the
    // parameters struct is copied per grid point/thread (see main_grid_d0.cpp),
    // so it cannot itself own a unique_ptr. The owning unique_ptr lives in
    // main() for the program's whole lifetime; only used for LHAPDF.
    Interpolator* D_frag_interp = nullptr;

    // Diffractive + fragmentation only: the gluon k is capped at the charm
    // quark's own momentum pc = pD0/zh -- the radiated gluon cannot carry
    // more transverse momentum than the charm quark that emitted it -- but
    // pc varies *within* the VEGAS call (zh is a sampled dimension), so that
    // cutoff can't be a VEGAS box edge directly. The fix: give k a VEGAS box
    // wide enough to cover every pc reachable at this pD0 -- k_upper_frag =
    // pD0/zmin, the largest pc possible over the whole zh range (since
    // pc = pD0/zh is largest at the smallest zh) -- then apply the exact
    // k<=pc cutoff as an explicit per-sample check inside the integrand
    // (same style as the x<=0||x>=1 kinematic checks elsewhere).
    double k_upper_frag;
    double log_k_frag;   // log(k_upper_frag / k_lo_frag), set before each VEGAS call
};

// D0-level (fragmented) exclusive integrand: dsigma/(d2pD0 dy), evaluated
// at the charm quark momentum pc=par->pD0/zh, with an extra fragmentation
// dimension zh in [par->zmin, par->zmax] and weight D(zh)/zh^2 folded in.
// 5D VEGAS: {r1, r2, u_qp, b, zh}.
double integrand_exclusive(double* vec, size_t dim, void* p);
double exclusiveCrossSection( void* p);

// D0-level (fragmented) diffractive integrand: dsigma/(d2pD0 dy), summed
// over x_po in [1e-6, 0.1] and over the gluon k, evaluated at pc=par->pD0/zh
// with the extra fragmentation dimension zh, and the gluon k capped at pc
// exactly (see k_upper_frag/log_k_frag above). 7D VEGAS:
// {r1, r2, u_k, u_qp, b, u_xpo, zh}.
double integrand_diffractive(double* vec, size_t dim, void* p);
double diffractiveCrossSection( void* p);

double flux_density(double qp, double b, void* p);
double photon_flux(double b, double qp, void* p);

#endif
