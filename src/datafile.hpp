/*
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

#ifndef _DATAFILE_H
#define _DATAFILE_H


#include <sstream>
#include <vector>

// BK solution read from a file
class DataFile
{
    public:
        DataFile(std::string fname);
        double MinR();
        double RMultiplier();
        int RPoints();
        
        double MaxY();
        double X0();

        // n[yind][rind] and the evolution rapidities
        void GetData(std::vector< std::vector<double> > &n,
            std::vector<double> &rapidities);

    private:
        std::string filename;
        std::vector<std::vector <double> > data;
        std::vector<double> yvals;
        double minr;
        double r_multiplier;
        int rpoints;
        double x0;

};

#endif
