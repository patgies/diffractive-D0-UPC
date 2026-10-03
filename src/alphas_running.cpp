#include "alphas_running.hpp"
#include <cmath>


static double alphas_mZ = 0.118;
static double mZ = 91.2;    
static double Nf = 4.0;     

double alphas_run(double mu)
{
    // one-loop running: 1/alpha_s(mu) = 1/alpha_s(mZ) + 2*b0*log(mu/mZ)
    double b0 = (33.0 - 2.0 * Nf) / (12.0 * M_PI);
    return 1.0 / (1.0 / alphas_mZ + 2.0 * b0 * log(mu / mZ));
}
