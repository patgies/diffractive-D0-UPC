#ifndef def_hpp
#define def_hpp
#include "amplitudelib.hpp"
#include "interpolation.hpp"
#include <memory>
#include <string>
#include <gsl/gsl_monte.h>

using namespace std;

enum class FragmentationType { BCFY, KniehlKramer, HymnD };

// Max x_po; dipole frozen above XPO_FROZEN (x0 of BK)
const double XPO_MAX = 0.1;
const double XPO_FROZEN = 0.01;
// Lower limits of the log sampling of k [GeV] and x_po
const double K_MIN = 0.01;
const double XPO_MIN = 1e-6;
// Max k when there is no k cut
const double K_MAX_NOCUT = 400.0;

struct parameters
{
    AmplitudeLib *dipole;
    double m;    // charm mass
    double y;
    double ss;
    // VEGAS calls. Exclusive: 5D below excl_pt_threshold, 3D above.
    size_t calls_excl;
    size_t calls_excl_factorized;
    double excl_pt_threshold;
    size_t calls_diff;
    double m2;   // m*m

    // Photon flux and nucleus
    double alpha, Z, mn, S;
    string channel;
    double bmin, bmax, qpmax;
    // "AA": Pb-Pb with Gamma_AA and EMD. "pA": proton target with Gamma_pA.
    string target = "AA";
    // "PL" point-like, "WS" Woods-Saxon, "EFF"/"STARLIGHT" effective (b already summed)
    string flux_model = "PL";

    // Fragmentation: pc = pD0/z_h, z_h in [z_h_min, z_h_max]
    double pD0;
    FragmentationType frag_type = FragmentationType::BCFY;
    double z_h_min, z_h_max;
    // D(z_h) at the fragmentation scale
    Interpolator* D_frag_interp = nullptr;

    // Diffractive cut k <= k_cut_factor * scale (K_CUT_FACTOR, <= 0: no cut).
    // Scale: pc for the quark DTMD; pc or m_T for the gluon DTMD (K_CUT_SCALE).
    double k_cut_factor = 2.0;
    string k_cut_scale = "mt";
    // Quark DTMD only for pc >= quark_p_min [GeV] (QUARK_P_MIN)
    double quark_p_min = 5.0;
    // "gluon", "quark" or "both" (DTMD)
    string dtmd = "both";
    double k_upper_frag;
    double log_k_frag;   // log(k_upper_frag / K_MIN)

    // Fixed x_po for the *_xpom functions
    double fixed_xpo;

    // Fixed-q+ mode: no flux, no fragmentation, p = charm pT
    double p, p2, mt;   // p, p^2, sqrt(p^2 + m^2)
    double kmax, k_upper, log_k;   // k range
};

// Effective fluxes have the b sum already done
inline bool is_effective_flux(const parameters* par)
{
    return par->flux_model == "EFF" || par->flux_model == "STARLIGHT";
}

// Lower limit of the b integral
inline double b_integration_min(const parameters* par)
{
    if (is_effective_flux(par)) return 0.0;
    if (par->target == "pA" || par->channel == "PL(AnAn)") return par->bmin;
    return 0.0;
}

// Settings from the environment (setup.cpp)
bool setup_parameters(parameters& param, unique_ptr<Interpolator>& d_frag_interp, double pD0, double y,
                      const string& frag_tag, const string& channel);

// VEGAS integral (vegas.cpp)
double vegas_integrate(gsl_monte_function* F, double* low, double* up, size_t calls);

// Exclusive radial integrals I_nu = int dr r J_nu(pc r) K_nu(m r) N(r), nu = 0, 1
void exclusive_radial_integrals(double pc, double m, double x_dip, AmplitudeLib* dipole, double& I0, double& I1);

// Exclusive dsigma/(d2pD0 dy): 5D VEGAS at low pt, 3D at high pt
double integrand_exclusive_mc(double* vec, size_t dim, void* p);
double integrand_exclusive_factorized(double* vec, size_t dim, void* p);
double exclusiveCrossSection( void* p);

// Exclusive at fixed x_po
double integrand_exclusive_xpom(double* vec, size_t dim, void* p);
double exclusiveCrossSection_xpom(void* p);

// Diffractive: 7D VEGAS
double integrand_diffractive(double* vec, size_t dim, void* p);
double diffractiveCrossSection( void* p);

// Diffractive at fixed x_po: 6D VEGAS
double integrand_diffractive_xpom(double* vec, size_t dim, void* p);
double diffractiveCrossSection_xpom(void* p);

// Flux (arXiv:2404.09731): r = photon to emitting nucleus, b = between nuclei,
// s = position in the target. b_d comes from the Glauber file, integrated in plotting_scripts/.
double flux_density(double qp, double r, void* p);
double flux_density_WS(double z_gamma, double r, void* p);
double photon_flux(double b, double qp, void* p);
// Effective flux (arXiv:2404.09731 Eq. 4)
double effective_photon_flux(double qp, void* p);

// Bessel functions
double Jn(int nu, double x);
double Kn(int nu, double x);

// Fixed q+ = 2p+, without prefactors (applied in charm_fixed_qp.py)
double exclusiveCrossSection_fixed_qp(void* p, double* q0 = nullptr, double* q1 = nullptr);
double integrand_diffractive_fixed_qp(double* vec, size_t dim, void* p);
double diffractiveCrossSection_fixed_qp(void* p);
// Same by quadrature, no Monte Carlo
double diffractiveCrossSection_fixed_qp_quad(void* p);
// int dr r J_nu(k r) K_nu(omega r) N(r) (N_A if adjoint), Gauss-Legendre
double dtmd_radial_integral(int nu, double k, double omega, double x_dip, bool adjoint, AmplitudeLib* dipole);
// Fixed x_po: dsigma/dln(x_po)
double integrand_diffractive_fixed_qp_xpom(double* vec, size_t dim, void* p);
double diffractiveCrossSection_fixed_qp_xpom(void* p);

#endif
