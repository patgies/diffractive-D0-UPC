#include "lhapdf_grid.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {

std::vector<double> ParseDoubles(const std::string& line)
{
    std::vector<double> out;
    std::istringstream iss(line);
    double v;
    while (iss >> v) out.push_back(v);
    return out;
}

std::vector<int> ParseInts(const std::string& line)
{
    std::vector<int> out;
    std::istringstream iss(line);
    int v;
    while (iss >> v) out.push_back(v);
    return out;
}

struct QGridBlock {
    std::vector<double> xgrid;
    std::vector<double> qgrid;
    std::vector<std::vector<double>> values;  // [ix][iq], already selected flavor column
};

}

std::unique_ptr<Interpolator> MakeLHAPDFZInterpolator(
    const std::string& member_file, int pdg_flavor, double Q)
{
    std::ifstream in(member_file);
    if (!in) throw std::runtime_error("MakeLHAPDFZInterpolator: cannot open " + member_file);

    std::string line;
    while (std::getline(in, line)) {
        if (line.compare(0, 3, "---") == 0) break;
    }

    std::vector<QGridBlock> blocks;

    while (std::getline(in, line)) {
        std::vector<double> xgrid = ParseDoubles(line);
        if (xgrid.empty()) break;

        if (!std::getline(in, line))
            throw std::runtime_error("MakeLHAPDFZInterpolator: truncated Q-grid line in " + member_file);
        std::vector<double> qgrid = ParseDoubles(line);

        if (!std::getline(in, line))
            throw std::runtime_error("MakeLHAPDFZInterpolator: truncated flavor line in " + member_file);
        std::vector<int> flavors = ParseInts(line);

        int flavor_col = -1;
        for (size_t i = 0; i < flavors.size(); i++) {
            if (flavors[i] == pdg_flavor) { flavor_col = (int)i; break; }
        }
        if (flavor_col < 0)
            throw std::runtime_error("MakeLHAPDFZInterpolator: flavor id not found in grid " + member_file);

        size_t nx = xgrid.size();
        size_t nq = qgrid.size();
        std::vector<std::vector<double>> values(nx, std::vector<double>(nq));

        for (size_t ix = 0; ix < nx; ix++) {
            for (size_t iq = 0; iq < nq; iq++) {
                if (!std::getline(in, line))
                    throw std::runtime_error("MakeLHAPDFZInterpolator: truncated value grid in " + member_file);
                std::vector<double> row = ParseDoubles(line);
                if ((int)row.size() <= flavor_col)
                    throw std::runtime_error("MakeLHAPDFZInterpolator: short value row in " + member_file);
                values[ix][iq] = row[flavor_col];
            }
        }

        std::getline(in, line);

        blocks.push_back({std::move(xgrid), std::move(qgrid), std::move(values)});
    }

    if (blocks.empty())
        throw std::runtime_error("MakeLHAPDFZInterpolator: no Q subgrids found in " + member_file);

    // LHAPDF subgrids are contiguous and increasing in Q; clamp an
    // out-of-range Q to the grid edges (frozen extrapolation, matching
    // standard PDF evolution behavior below Q0/above Qmax) instead of
    // failing hard. This matters for scale-variation studies (e.g.
    // Q = mt0/2) that can dip below QMin for low-pT points.
    double q_clamped = Q;
    if (q_clamped < blocks.front().qgrid.front()) q_clamped = blocks.front().qgrid.front();
    if (q_clamped > blocks.back().qgrid.back()) q_clamped = blocks.back().qgrid.back();

    for (const auto& block : blocks) {
        const std::vector<double>& qgrid = block.qgrid;
        if (q_clamped < qgrid.front() || q_clamped > qgrid.back()) continue;

        size_t nx = block.xgrid.size();
        size_t iq_hi = 1;
        while (iq_hi < qgrid.size() - 1 && qgrid[iq_hi] < q_clamped) iq_hi++;
        size_t iq_lo = iq_hi - 1;
        double q_lo = qgrid[iq_lo], q_hi = qgrid[iq_hi];
        double t = (q_hi > q_lo) ? (q_clamped - q_lo) / (q_hi - q_lo) : 0.0;

        std::vector<double> zvals(nx), dvals(nx);
        for (size_t ix = 0; ix < nx; ix++) {
            double xf = block.values[ix][iq_lo] * (1.0 - t) + block.values[ix][iq_hi] * t;
            zvals[ix] = block.xgrid[ix];
            dvals[ix] = xf / zvals[ix];   // grid stores x*f(x,Q) -> convert to D(z,Q)
        }

        return std::unique_ptr<Interpolator>(new Interpolator(zvals, dvals));
    }

    throw std::runtime_error("MakeLHAPDFZInterpolator: clamped Q not found in any subgrid in " + member_file);
}
