import glob
import math
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson
sys.path.insert(0, os.path.dirname(__file__))
from alphas_running import alphas_run

# Integrates the per-b, per-(pD0,y) grid files written by
# ../local_workflows/run_many_nucleus.sh (src/main.cpp / D0) over the target-nucleus
# impact parameter b (Simpson's rule, weighted by b)
#   b column = target-nucleus Glauber sample (data/<NUCLEUS>/mve/glauber_mve_<b>)
#   exclusive:   p_e = alpha_em * Nc * e_c^2 / (2*pi^2)
#   diffractive: p_d = alpha_s(mT) * alpha_em * e_c^2 * (Nc^2-1) / (8*pi^4)
# (no sigma0 -- see prefactor()'s docstring) and 2*pi*pD0 Jacobian

alphae = 1/137
mc     = 1.5          # charm mass in GeV
e_c    = 2/3
Nc     = 3

# Unit conversion GeV^-2 -> mb 
FMGEV = 5.068
GEVSQR_TO_NB = 1.0e7 / (FMGEV * FMGEV)
GEVSQR_TO_MB = GEVSQR_TO_NB * 1e-6

NUCLEUS = os.environ.get("NUCLEUS", "Pb")
FRAG    = os.environ.get("FRAG_TYPE", "BCFY")
CHANNEL = os.environ.get("CHANNEL", "An0n").translate(str.maketrans('', '', '() '))


def read_rapidity(filename):
    with open(filename) as f:
        for line in f:
            if "fixed rapidity y" in line:
                return float(line.split(":")[-1])
    raise ValueError(f"no rapidity header found in {filename}")


def read_data_file(filename):
    """Returns three plain lists (b, pD0, dsigma), skipping comment lines."""
    b_list, pt_list, dsigma_list = [], [], []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            b, pt, dsigma = line.split()
            b_list.append(float(b))
            pt_list.append(float(pt))
            dsigma_list.append(float(dsigma))
    return b_list, pt_list, dsigma_list


def group_by_pt(b_list, pt_list, dsigma_list):
    groups = {}
    for b, pt, dsigma in zip(b_list, pt_list, dsigma_list):
        if pt not in groups:
            groups[pt] = []
        groups[pt].append((b, dsigma))
    return groups


def integrate_over_b(pairs):
    """int b*dsigma(b) db via Simpson's rule (the radial Glauber-b measure)."""
    pairs = sorted(pairs)   # sorts by b first since these are (b, dsigma) tuples
    if pairs[0][0] < 0:
        raise ValueError(
            "Negative b found -- datafile wasn't a Glauber sample "
            "(e.g. the plain proton dipole has no b to integrate over)."
        )
    b_values = [pair[0] for pair in pairs]
    weighted = [b * dsigma for b, dsigma in pairs]
    return simpson(weighted, x=b_values)


def prefactor(process, pt):
    """Physical prefactor, matching D0.py's convention.
    No sigma0 here: that's the GBW proton-dipole normalization, and for a
    nucleus target the "area" is already accounted for by integrate_over_b's
    Glauber b-integral -- applying sigma0 on top of that would double-count it.
    """
    if process == "exclusive":
        return alphae * Nc * e_c**2 / (2 * math.pi**2)
    elif process == "diffractive":
        alphas = alphas_run(math.sqrt(pt**2 + mc**2))
        return alphas * alphae * e_c**2 * (Nc**2 - 1) / (8 * math.pi**4)
    raise ValueError(f"unknown process {process}")


def load_results(process):
    pattern = f"../output/D0_{process}_{FRAG}_{CHANNEL}_{NUCLEUS}_y*.dat"
    results = {}
    for filename in sorted(glob.glob(pattern)):
        y = read_rapidity(filename)
        b_list, pt_list, dsigma_list = read_data_file(filename)
        pt_groups = group_by_pt(b_list, pt_list, dsigma_list)

        results[y] = []
        for pt, pairs in sorted(pt_groups.items()):
            b_integral = integrate_over_b(pairs)
            # 2*pi*b_integral: Glauber transverse-plane (b) integral.
            # 2*pi*pt: d^2pD0 -> dpD0 "kt-spectrum" Jacobian (see plot_pt_spectrum.py).
            cross_section = (2*math.pi) * b_integral * prefactor(process, pt) \
                            * (2*math.pi) * pt * GEVSQR_TO_MB
            results[y].append((pt, cross_section))
    return results


def write_combined(process, results):
    for y, points in sorted(results.items()):
        ytag = str(y).replace(".", "")
        outfile = f"../output/combined_d0_{process}_{FRAG}_{CHANNEL}_{NUCLEUS}_y{ytag}.dat"
        with open(outfile, "w") as f:
            f.write(f"# {process} D0 cross section, {NUCLEUS} target, {FRAG} fragmentation, "
                    f"{CHANNEL} channel, integrated over Glauber b\n")
            f.write(f"# fixed rapidity y : {y}\n")
            f.write("# pD0  dsigma_dy_dpD0 [mb/GeV]\n")
            for pt, cs in points:
                f.write(f"{pt:.6g}  {cs:.6g}\n")
        print(f"Wrote {outfile}")


def main():
    results_excl = load_results("exclusive")
    results_diff = load_results("diffractive")
    if not results_excl and not results_diff:
        sys.exit(f"No files found for NUCLEUS={NUCLEUS}, FRAG_TYPE={FRAG}, CHANNEL={CHANNEL} "
                  "in ../output/ -- run ../local_workflows/run_many_nucleus.sh first.")

    write_combined("exclusive", results_excl)
    write_combined("diffractive", results_diff)

    y_values = sorted(set(results_excl) | set(results_diff))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(y_values)))

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for y, c in zip(y_values, colors):
        if y in results_diff:
            points = sorted(results_diff[y])
            pts = [point[0] for point in points]
            cs = [point[1] for point in points]
            ax.plot(pts, cs, color=c, linestyle='-', label=f"diff., $y={y:g}$")
        if y in results_excl:
            points = sorted(results_excl[y])
            pts = [point[0] for point in points]
            cs = [point[1] for point in points]
            ax.plot(pts, cs, color=c, linestyle='--', label=f"excl., $y={y:g}$")

    ax.set_yscale('log')
    ax.set_xlabel(r'$p_{D^0\perp}$ [GeV]')
    ax.set_ylabel(r'$d\sigma/dy\,dp_{D^0\perp}$ (mb/GeV)')
    ax.set_title(f'D0-level $p_{{D^0}}$ spectrum, {NUCLEUS} target ({FRAG}, {CHANNEL})')
    ax.tick_params(axis='both', which='major', labelsize=12)

    handles, labels = ax.get_legend_handles_labels()
    diff_items = [(h, l) for h, l in zip(handles, labels) if l.startswith('diff')]
    excl_items = [(h, l) for h, l in zip(handles, labels) if l.startswith('excl')]
    handles, labels = zip(*(diff_items + excl_items))
    ax.legend(handles, labels, fontsize=7.5, ncol=2)
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    outname = f"../plots/combined_d0_{FRAG}_{CHANNEL}_{NUCLEUS}{NUCLEUS}.png"
    plt.savefig(outname, dpi=300)
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
