import glob
import math
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.integrate import simpson
sys.path.insert(0, os.path.dirname(__file__))
from alphas_running import alphas_run

# This script reads the data files written by ../run_many_nucleus.sh
# (produced by the D0 program), integrates them over the target-nucleus
# impact parameter b using Simpson's rule, and makes one plot with every
# combination of process (exclusive/diffractive) and fragmentation function
# (BCFY/KniehlKramer) on it. Each rapidity gets its own color, and each
# (process, fragmentation) combination gets its own linestyle.

alphae = 1/137
mc     = 1.5          # charm mass in GeV
e_c    = 2/3
Nc     = 3
sigma0 = 16.36         # mb

# conversion factor from GeV^-2 to mb
FMGEV = 5.068
GEVSQR_TO_NB = 1.0e7 / (FMGEV * FMGEV)
GEVSQR_TO_MB = GEVSQR_TO_NB * 1e-6

NUCLEUS = os.environ.get("NUCLEUS", "Pb")
CHANNEL = os.environ.get("CHANNEL", "An0n").translate(str.maketrans('', '', '() '))

PROCESSES = ["diffractive", "exclusive"]
FRAG_TYPES = ["BCFY", "KniehlKramer"]

# which linestyle to use for each (process, fragmentation) combination
LINESTYLES = {
    ("diffractive", "BCFY"):          "-",
    ("diffractive", "KniehlKramer"):  "--",
    ("exclusive",   "BCFY"):          ":",
    ("exclusive",   "KniehlKramer"):  "-.",
}


def read_rapidity(filename):
    with open(filename) as f:
        for line in f:
            if "fixed rapidity y" in line:
                return float(line.split(":")[-1])
    raise ValueError(f"no rapidity header found in {filename}")


def read_data_file(filename):
    """Returns three plain lists (b, pD0, dsigma), skipping comment lines."""
    b_list = []
    pt_list = []
    dsigma_list = []
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
    """Physical prefactor, matching plot_pt_spectrum.py's convention."""
    if process == "exclusive":
        return alphae * Nc * e_c**2 * sigma0 / (2 * math.pi**2)
    elif process == "diffractive":
        alphas = alphas_run(math.sqrt(pt**2 + mc**2))
        return alphas * alphae * e_c**2 * (Nc**2 - 1) * sigma0 / (8 * math.pi**4)
    raise ValueError(f"unknown process {process}")


def load_results(process, frag):
    """Read all files for one (process, frag) combination and return {y: [(pt, cross_section), ...]}."""
    pattern = f"../files/d0_point_{process}_{frag}_{CHANNEL}_{NUCLEUS}_y*.dat"
    results = {}
    for filename in sorted(glob.glob(pattern)):
        y = read_rapidity(filename)
        b_list, pt_list, dsigma_list = read_data_file(filename)
        pt_groups = group_by_pt(b_list, pt_list, dsigma_list)

        results[y] = []
        for pt, pairs in sorted(pt_groups.items()):
            b_integral = integrate_over_b(pairs)
            # 2*pi*b_integral is the Glauber transverse-plane (b) integral.
            # 2*pi*pt is the Jacobian from d^2pD0 to dpD0.
            cross_section = (2*math.pi) * b_integral * prefactor(process, pt) \
                            * (2*math.pi) * pt * GEVSQR_TO_MB
            results[y].append((pt, cross_section))
    return results


def main():
    all_results = {}
    for process in PROCESSES:
        for frag in FRAG_TYPES:
            results = load_results(process, frag)
            if results:
                all_results[(process, frag)] = results

    if not all_results:
        sys.exit(f"No files found for NUCLEUS={NUCLEUS}, CHANNEL={CHANNEL} "
                  "in ../files/ -- run ../run_many_nucleus.sh first.")

    # collect every rapidity value that shows up in any of the results,
    # keeping only y <= 2.0 so the plot doesn't get too crowded
    y_values = set()
    for results in all_results.values():
        for y in results:
            if y <= 2.0:
                y_values.add(y)
    y_values = sorted(y_values)

    color_list = plt.cm.viridis(np.linspace(0.15, 0.85, len(y_values)))
    colors = {}
    for y, c in zip(y_values, color_list):
        colors[y] = c

    plt.figure(figsize=(7.5, 6.5))

    for process, frag in LINESTYLES:
        results = all_results.get((process, frag))
        if not results:
            continue
        linestyle = LINESTYLES[(process, frag)]
        for y in sorted(results):
            if y > 2.0:
                continue
            points = sorted(results[y])
            pt_values = [pair[0] for pair in points]
            cross_section_values = [pair[1] for pair in points]
            plt.plot(pt_values, cross_section_values, color=colors[y], linestyle=linestyle)

    plt.yscale("log")
    plt.xlabel(r"$p_{D^0}$ [GeV]")
    plt.ylabel(r"$d\sigma/dy\,dp_{D^0}$ [mb/GeV]")
    plt.title(f"$D^0$ photoproduction, {NUCLEUS}+{NUCLEUS} UPC ({CHANNEL})")

    y_handles = [Line2D([0], [0], color=colors[y], linestyle="-", label=f"y={y:g}") for y in y_values]
    y_legend = plt.legend(handles=y_handles, loc="upper right")
    plt.gca().add_artist(y_legend)

    style_handles = []
    for (process, frag), linestyle in LINESTYLES.items():
        if (process, frag) in all_results:
            style_handles.append(Line2D([0], [0], color="black", linestyle=linestyle, label=f"{process}, {frag}"))
    plt.legend(handles=style_handles, loc="lower left", fontsize=8)

    plt.tight_layout()
    outname = f"../plots/cross_section_{CHANNEL}_{NUCLEUS}.png"
    plt.savefig(outname, dpi=150)
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
