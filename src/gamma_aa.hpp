#ifndef GAMMA_AA_HPP
#define GAMMA_AA_HPP
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <cmath>
#include <string>
#include <vector>
#include <stdexcept>
#include <fstream>
 
 


class TASpline {
private:
    gsl_spline* spline;
    gsl_interp_accel* acc;
    double s_min, s_max;
 
public:
    TASpline(const std::vector<double>& s,
             const std::vector<double>& T)
    {
        spline = gsl_spline_alloc(gsl_interp_cspline, s.size());
        acc = gsl_interp_accel_alloc();
        gsl_spline_init(spline, s.data(), T.data(), s.size());
 
        s_min = s.front();
        s_max = s.back();
    }
 
    
    TASpline(const TASpline&) = delete;
    TASpline& operator=(const TASpline&) = delete;
 
    double operator()(double s)
    {
        if (s < s_min || s > s_max)
            return 1.0;
 
        return gsl_spline_eval(spline, s, acc);
    }
 
    ~TASpline()
    {
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
    }
};
 

class GammaAA {
private:
    TASpline& TA;
    double b_min;
    double b_max;
 
public:
    GammaAA(TASpline& ta)
        : TA(ta), b_min(0.0), b_max(1e6)
    {}
 
    double operator()(double b)
    {
        if (b < b_min)
            throw std::out_of_range("b too small");
 
        if (b > b_max)
            return 1.0;
 
        return TA(b); // Directly calls the spline
    }
};
 
 

void load_data_and_initialize(const std::string& filename);
 

GammaAA& gamma_aa();
 
#endif
