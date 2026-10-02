#ifndef KK_GRID_HPP
#define KK_GRID_HPP

#include <memory>
#include "interpolation.hpp"

// Kniehl-Kramer c -> D0 fragmentation function D(z), evolved with DGLAP,
// at the scale Q (GeV).

std::unique_ptr<Interpolator> MakeKniehlKramerInterpolator(double Q);

#endif
