#include "kk_grid.hpp"
#include "hymnd_grid.hpp"
#include <cstdlib>

namespace {
const char* DefaultGridPath() { return "inputs/kk_eko/kk_eko_0000.dat"; }
const int kCharmPid = 4;
}

std::unique_ptr<Interpolator> MakeKniehlKramerInterpolator(double Q)
{
    // KK_EKO_FILE lets a member/replica grid override the default, same
    // convention as HYMND_FILE for frag_type=HymnD.
    const char* path = std::getenv("KK_EKO_FILE");
    return MakeHymnDZInterpolator(path ? path : DefaultGridPath(), kCharmPid, Q);
}
