#include "kk_grid.hpp"
#include "hymnd_grid.hpp"
#include <cstdlib>

namespace {
const char* DefaultGridPath() { return "inputs/KK_EKO/kk_eko_0000.dat"; }
const int kCharmPid = 4;
}

std::unique_ptr<Interpolator> MakeKniehlKramerInterpolator(double Q)
{
    // KK_EKO_FILE: use another grid file.
    const char* path = std::getenv("KK_EKO_FILE");
    return MakeHymnDZInterpolator(path ? path : DefaultGridPath(), kCharmPid, Q);
}
