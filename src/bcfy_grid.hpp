#ifndef BCFY_GRID_HPP
#define BCFY_GRID_HPP

#include <memory>
#include "interpolation.hpp"

// BCFY c -> D0 fragmentation function D(z_h), evolved with DGLAP,
// at the scale Q (GeV).

std::unique_ptr<Interpolator> MakeBCFYInterpolator(double Q);

#endif
