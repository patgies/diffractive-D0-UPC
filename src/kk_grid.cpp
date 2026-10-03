#include "kk_grid.hpp"
#include "lhapdf_grid.hpp"
#include <cstdlib>

namespace {
const char* DefaultGridPath() { return "input/KK_EKO/kk_eko_0000.dat"; }
const int kCharmPid = 4;
}

std::unique_ptr<Interpolator> MakeKniehlKramerInterpolator(double Q)
{
    const char* path = std::getenv("KK_EKO_FILE");
    return MakeLHAPDFGridInterpolator(path ? path : DefaultGridPath(), kCharmPid, Q);
}
