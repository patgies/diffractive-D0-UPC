#ifndef BCFY_GRID_HPP
#define BCFY_GRID_HPP

#include <memory>
#include "interpolation.hpp"

// BCFY c -> D0 fragmentation D(z_h, Q), DGLAP-evolved

std::unique_ptr<Interpolator> MakeBCFYInterpolator(double Q);

#endif
