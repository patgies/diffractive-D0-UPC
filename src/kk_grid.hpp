#ifndef KK_GRID_HPP
#define KK_GRID_HPP

#include <memory>
#include "interpolation.hpp"

// Kniehl-Kramer c -> D0 fragmentation function D(z_h), evolved with DGLAP

std::unique_ptr<Interpolator> MakeKniehlKramerInterpolator(double Q);

#endif
