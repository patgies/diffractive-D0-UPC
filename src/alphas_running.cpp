#include "alphas_running.hpp"
#include <cmath>

// PDG value of the strong coupling at the Z boson mass
static double alphas_mZ = 0.118;
static double mZ = 91.2;    // Z boson mass, in GeV
static double Nf = 4.0;     // number of active quark flavors at the charm mass scale

double alphas_run(double mu)
{
    // One-loop running coupling formula:
    //   1/alpha_s(mu) = 1/alpha_s(mZ) + 2*b0*log(mu/mZ)
    double b0 = (33.0 - 2.0 * Nf) / (12.0 * M_PI);
    return 1.0 / (1.0 / alphas_mZ + 2.0 * b0 * log(mu / mZ));
}
