#ifndef KK_GRID_HPP
#define KK_GRID_HPP

#include <memory>
#include "interpolation.hpp"

// Builds a z-interpolator for the DGLAP-evolved Kniehl & Kramer c -> D0
// fragmentation function at a fixed factorisation scale Q (GeV).

std::unique_ptr<Interpolator> MakeKniehlKramerInterpolator(double Q);

#endif
