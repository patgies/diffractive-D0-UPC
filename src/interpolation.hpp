/*
 * General purpose interpolation class
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

#ifndef _INTERPOLATION_H
#define _INTERPOLATION_H



#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <vector>
#include <iostream>

using std::cout;
using std::endl;
using std::cerr;


enum INTERPOLATION_METHOD {
    INTERPOLATE_SPLINE
};


// GSL spline interpolation
class Interpolator
{
    public:
        // From two arrays of p points (pointers are kept)
        Interpolator(double* x, double* y, int p);

        // From two vectors (copied)
        Interpolator(std::vector<double> &x, std::vector<double> &y);
        Interpolator(const Interpolator& inter);
        ~Interpolator();
        void Clear();
        double Evaluate(double x);
        // f'(x)
        double Derivative(double x);    
        // f''(x)
        double Derivative2(double x);   
        void SetMethod(INTERPOLATION_METHOD m);

        // Called by the constructor; call again after SetMethod()
        int Initialize();

        double MinX();
        double MaxX();

        double* GetXData();
        double* GetYData();
        gsl_spline* GetGslSpline() const;
        int GetNumOfPoints() const;
        INTERPOLATION_METHOD GetMethod() const;
        
        bool Freeze();

        // Outside the range return the under/overflow values
        void SetFreeze(bool f);
        void SetUnderflow(double min);
        void SetOverflow(double max);
        double UnderFlow();
        double OverFlow();

        // Interpolation range (default: the data range)
        void SetMaxX(double x);
        void SetMinX(double x);

        // Print an error when evaluated outside the range
        void SetOutOfRangeErrors(bool set);

    private:
        INTERPOLATION_METHOD method;
        double* xdata, *ydata;
        bool allocated_data;    // we own xdata, ydata
        double minx,maxx;
        int points;
        bool ready;
        
        bool freeze;		// return the under/overflow values outside the range
        double freeze_underflow;
		double freeze_overflow;
        
        gsl_interp_accel *acc;
        gsl_spline *spline;

      
        bool out_of_range_errors;


};




#endif
