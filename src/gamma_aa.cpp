#include "gamma_aa.hpp"
#include "def.hpp"
#include <memory>
#include <cmath>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

static double Kn_flux(int nu, double x)
{
    gsl_sf_result result;
    int status = gsl_sf_bessel_Kn_e(nu, x, &result);
    if (status != GSL_SUCCESS) return 0.0;
    return result.val;
}

static double get_GammaAA(double b) { return gamma_aa()(b); }

double flux_density(double qp, double b, void* p)
{
    parameters* par = (parameters*)p;
    double omega = qp / sqrt(2.0);
    double z     = omega * 2.0 / par->ss;
    double eta   = z * par->mn * b;
    if (eta > 50.0) return 0.0;
    double pref = (par->alpha * par->Z * par->Z) / (M_PI * M_PI);
    double K1   = Kn_flux(1, eta);
    return (pref / omega) * (eta*eta / (b*b)) * (K1*K1);
}

double photon_flux(double b, double qp, void* p)
{
    parameters* par = (parameters*)p;
    const string& channel = par->channel;
    double flux  = flux_density(qp, b, p);
    double gamma = (channel == "PL(AnAn)") ? ((b >= 14.0/0.197327) ? 1.0 : 0.0)
                                           : ((b < 150.0) ? get_GammaAA(b) : 1.0);
    double P_emd      = par->S / (b * b);
    double P_no_emd   = std::exp(-P_emd);
    double emd_factor = 1.0;
    if      (channel == "An0n") emd_factor = P_no_emd;
    else if (channel == "Xn0n") emd_factor = P_no_emd * (1.0 - P_no_emd);
    return 2.0 * M_PI * b * flux * gamma * emd_factor;
}


namespace {
    std::unique_ptr<TASpline> ta_spline;
    std::unique_ptr<GammaAA>  gamma_instance;
}

void load_data_and_initialize(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file!");
    }

    std::vector<double> b_values;
    std::vector<double> T_values;
    double b, T;

    while (file >> b >> T) {
        b_values.push_back(b);
        T_values.push_back(T);
    }


    ta_spline.reset(new TASpline(b_values, T_values));
    gamma_instance.reset(new GammaAA(*ta_spline));
}

GammaAA& gamma_aa()
{
    if (!gamma_instance) {
        throw std::runtime_error(
            "gamma_aa() called before load_data_and_initialize()");
    }
    return *gamma_instance;
}