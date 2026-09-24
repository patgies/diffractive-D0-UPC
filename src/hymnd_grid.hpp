#ifndef HYMND_GRID_HPP
#define HYMND_GRID_HPP

#include <memory>
#include <string>
#include "interpolation.hpp"


std::unique_ptr<Interpolator> MakeHymnDZInterpolator(
    const std::string& member_file, int pdg_flavor, double Q);

#endif
