#ifndef BCFY_GRID_HPP
#define BCFY_GRID_HPP

#include <memory>
#include "interpolation.hpp"

// Builds a z-interpolator for the DGLAP-evolved Braaten-Cheung-Fleming-Yuan
// (BCFY) c -> D0 fragmentation function at a fixed factorisation scale Q
// (GeV). 

std::unique_ptr<Interpolator> MakeBCFYInterpolator(double Q);

#endif
