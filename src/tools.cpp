/*
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2015
 */

#include "tools.hpp"
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <gsl/gsl_integration.h>



double StrToReal(std::string str)
{
    std::stringstream buff(str);
    double tmp;
    buff >> tmp;
    return tmp;
}

int StrToInt(std::string str)
{
    std::stringstream buff(str);
    int tmp;
    buff >> tmp;
    return tmp;
}

// Index i with vec[i] <= val < vec[i+1] (vec increasing); -1 below the range

int FindIndex(double val, std::vector<double> &vec)
{
    if (val < vec[0]) return -1;
    
    int ind=-1;
    
    uint start=0; uint end=vec.size()-1;
    while(end-start>5)
    {
        int tmp = static_cast<int>((start+end)/2.0);
        
        if (vec[tmp]>=val)
            end=tmp;
        else
            start=tmp;
    }
    
    
    for (uint i=start; i<=end; i++)
    {
        if (vec[i]<=val and vec[i+1]>val)
        {
            ind=i;
            break;
        }
    }
    if (ind == -1) return vec.size()-1;
    return ind;
}

