#ifndef KK_GRID_HPP
#define KK_GRID_HPP

#include <memory>
#include "interpolation.hpp"

// Kniehl-Kramer c -> D0 fragmentation D(z_h, Q), DGLAP-evolved

std::unique_ptr<Interpolator> MakeKniehlKramerInterpolator(double Q);

#endif
