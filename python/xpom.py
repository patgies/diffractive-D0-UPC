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


plt.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{amssymb}",  # for \mathbb{P} below
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 28,
    "axes.titlesize": 20,
    "xtick.labelsize": 24,
    "ytick.labelsize": 24,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 8,
    "ytick.major.size": 8,
    "xtick.minor.size": 4,
    "ytick.minor.size": 4,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    # same gray ticks with black numbers as D0.py
    "xtick.color": "0.4",
    "ytick.color": "0.4",
    "xtick.labelcolor": "black",
    "ytick.labelcolor": "black",
    "xtick.major.pad": 8,
    "ytick.major.pad": 4,
})

# This script reads the data files written by ../local_workflows/run_many_xpom.sh (produced
# by the D0_xpom program), integrates them over the target-nucleus impact
# parameter b using Simpson's rule, and plots the diffractive cross section
# as a function of ln(x_po), for each rapidity and pD0 value.

alphae = 1/137
mc     = 1.5          # charm mass in GeV
e_c    = 2/3
Nc     = 3

# conversion factor from GeV^-2 to mb
FMGEV = 5.068
GEVSQR_TO_NB = 1.0e7 / (FMGEV * FMGEV)
GEVSQR_TO_MB = GEVSQR_TO_NB * 1e-6

NUCLEUS = os.environ.get("NUCLEUS", "Pb")
FRAG    = os.environ.get("FRAG_TYPE", "BCFY")
CHANNEL = os.environ.get("CHANNEL", "An0n").translate(str.maketrans('', '', '() '))


def read_header(filename):
    """Returns (y, pD0) read from the '# fixed rapidity y :' / '# fixed pD0 :' header lines."""
    y = None
    pt = None
    with open(filename) as f:
        for line in f:
            if "fixed rapidity y" in line:
                y = float(line.split(":")[-1])
            elif "fixed pD0" in line:
                pt = float(line.split(":")[-1])
    if y is None or pt is None:
        raise ValueError(f"missing rapidity/pD0 header in {filename}")
    return y, pt


def read_data_file(filename):
    """Returns three plain lists (b, x_po, dsigma), skipping comment lines."""
    b_list = []
    xpo_list = []
    dsigma_list = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            b, xpo, dsigma = line.split()
            b_list.append(float(b))
            xpo_list.append(float(xpo))
            dsigma_list.append(float(dsigma))
    return b_list, xpo_list, dsigma_list


def group_by_xpo(b_list, xpo_list, dsigma_list):
    groups = {}
    for b, xpo, dsigma in zip(b_list, xpo_list, dsigma_list):
        if xpo not in groups:
            groups[xpo] = []
        groups[xpo].append((b, dsigma))
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


def diffractive_prefactor(pt):
    """Physical prefactor, matching D0.py's convention. No sigma0:
    that's the GBW proton-dipole normalization, and for a nucleus target the
    "area" is already accounted for by integrate_over_b's Glauber b-integral
    -- applying sigma0 on top of that would double-count it."""
    alphas = alphas_run(math.sqrt(pt**2 + mc**2))
    return alphas * alphae * e_c**2 * (Nc**2 - 1) / (8 * math.pi**4)


def exclusive_prefactor():
    """Physical prefactor for exclusive, matching diffractive_prefactor above
    -- no pt/alpha_s dependence for exclusive, and no sigma0 for the same
    reason."""
    return alphae * Nc * e_c**2 / (2 * math.pi**2)


def load_results(process="diffractive"):
    """Read all files for one process ("exclusive"/"diffractive") and return
    {(y, pt): [(x_po, dsigma_dy_dpt_dlnxpo), ...]}."""
    prefactor = diffractive_prefactor if process == "diffractive" else (lambda pt: exclusive_prefactor())
    pattern = f"../output/{CHANNEL}/xpom/D0_{process}_xpom_{FRAG}_{CHANNEL}_{NUCLEUS}_y*_pt*.dat"
    results = {}
    for filename in sorted(glob.glob(pattern)):
        y, pt = read_header(filename)
        b_list, xpo_list, dsigma_list = read_data_file(filename)
        xpo_groups = group_by_xpo(b_list, xpo_list, dsigma_list)

        points = []
        for xpo, pairs in sorted(xpo_groups.items()):
            b_integral = integrate_over_b(pairs)
            # 2*pi*b_integral is the Glauber transverse-plane (b) integral.
            # 2*pi*pt is the Jacobian from d^2pD0 to dpD0.
            dsigma_dxpo = (2*math.pi) * b_integral * prefactor(pt) \
                          * (2*math.pi) * pt * GEVSQR_TO_MB
            # dln(x_po) = dx_po / x_po
            dsigma_dlnxpo = xpo * dsigma_dxpo
            points.append((xpo, dsigma_dlnxpo))
        results[(y, pt)] = sorted(points)
    return results


def main():
    results = load_results("diffractive")
    if not results:
        sys.exit(f"No files found for NUCLEUS={NUCLEUS}, FRAG_TYPE={FRAG}, CHANNEL={CHANNEL} "
                  f"in ../output/{CHANNEL}/xpom/ -- run ../local_workflows/run_many_xpom.sh first.")

    # pD0 values to include in the plot:
    PT_TO_PLOT = [
        # 1.0,
        2.0,
        # 4.0,
        # 8.0,
    ]
    filtered_results = {}
    for (y, pt), points in results.items():
        if pt in PT_TO_PLOT:
            filtered_results[(y, pt)] = points
    results = filtered_results

    y_values = sorted(set(y for y, pt in results))
    pt_values = sorted(set(pt for y, pt in results))

    # Blue -> red steps (ColorBrewer RdBu 8-class without its two palest
    # middle steps), one per rapidity in increasing y -- the 6-y extension of
    # D0.py's 4-step palette. Linestyle also changes with y, so
    # curves stay distinguishable where they overlap and in black-and-white.
    palette = ["#2166ac", "#4393c3", "#92c5de", "#f4a582", "#d6604d", "#b2182b"]
    if len(y_values) <= 4:
        palette = ["#2166ac", "#67a9cf", "#ef8a62", "#b2182b"]   # D0.py's
    colors = {y: palette[round(i * (len(palette) - 1) / max(len(y_values) - 1, 1))]
              for i, y in enumerate(y_values)}
    linestyle_cycle = ["-", "--", ":", "-.", (0, (5, 1.5, 1, 1.5, 1, 1.5)), (0, (8, 2))]
    y_linestyles = {y: linestyle_cycle[i % len(linestyle_cycle)] for i, y in enumerate(y_values)}

    plt.figure(figsize=(8, 7))
    for (y, pt), points in sorted(results.items()):
        xpo_values = [p[0] for p in points]
        dsigma_values = [p[1] for p in points]
        ln_xpo = np.log(xpo_values)
        plt.plot(ln_xpo, dsigma_values, color=colors[y], linestyle=y_linestyles[y], lw=2.0)

    plt.yscale("log")
    plt.ylim(1e-9, 1e2)   # headroom above the curves (max ~0.2) for the legend
    plt.xlim(left=-12)
    plt.xlabel(r"$\ln(x_{\mathbb{P}})$", labelpad=14)
    plt.ylabel(r"$d\sigma/dy\,dp_{D^0\perp}\,d\ln(x_{\mathbb{P}})$ [mb/GeV]", labelpad=16)

    # Secondary top axis showing x_po 
    ax = plt.gca()
    ax.tick_params(axis='x', which='both', top=False)  # primary axis's own top ticks would clash with secax's
    secax = ax.secondary_xaxis('top', functions=(np.exp, np.log))
    secax.set_xscale('log')
    secax.set_xlabel(r"$x_{\mathbb{P}}$", labelpad=10)

    # in-plot label instead of a title, as in D0.py
    process_label = "Diffractive"
    pt_label = ", ".join(f"{pt:g}" for pt in pt_values)
    # bottom right: the spectra rise toward large x_po, so the top right is taken
    ax.text(0.95, 0.05,
            f"{NUCLEUS}-{NUCLEUS} 5.36 TeV\n{CHANNEL}, {FRAG}\n{process_label}\n"
            f"$p_{{D^0\\perp}} = {pt_label}$ GeV",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=20, linespacing=1.4)

    y_handles = [Line2D([0], [0], color=colors[y], linestyle=y_linestyles[y], lw=2.0, label=f"$y={y:g}$")
                 for y in y_values]
    plt.legend(handles=y_handles, loc="upper left", ncol=2, fontsize=20, frameon=False, columnspacing=1.2)

    plt.tight_layout()
    outname = f"../plots/xpom_{CHANNEL}_{NUCLEUS}{NUCLEUS}.pdf"
    plt.savefig(outname, dpi=150)
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
