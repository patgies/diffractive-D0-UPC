#ifndef LHAPDF_GRID_HPP
#define LHAPDF_GRID_HPP

#include <memory>
#include <string>
#include "interpolation.hpp"


std::unique_ptr<Interpolator> MakeLHAPDFZInterpolator(
    const std::string& member_file, int pdg_flavor, double Q);

#endif
