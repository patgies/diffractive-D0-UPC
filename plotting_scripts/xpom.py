import glob
import math
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

# Reads the D0_xpom output of run_xpom.sh, integrates over b_d and plots
# the diffractive cross section versus ln(x_po).
from D0 import integrate_over_b_d, prefactor, GEVSQR_TO_MB, NUCLEUS, CHANNEL   # also sets the plot style

plt.rcParams["text.latex.preamble"] = r"\usepackage{amssymb}"   # \mathbb{P}

FRAG = os.environ.get("FRAG_TYPE", "BCFY")


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
    """Returns three plain lists (b_d, x_po, dsigma), skipping comment lines."""
    b_d_list = []
    xpo_list = []
    dsigma_list = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            b_d, xpo, dsigma = line.split()
            b_d_list.append(float(b_d))
            xpo_list.append(float(xpo))
            dsigma_list.append(float(dsigma))
    return b_d_list, xpo_list, dsigma_list


def group_by_xpo(b_d_list, xpo_list, dsigma_list):
    groups = {}
    for b_d, xpo, dsigma in zip(b_d_list, xpo_list, dsigma_list):
        if xpo not in groups:
            groups[xpo] = []
        groups[xpo].append((b_d, dsigma))
    return groups


def load_results(process="diffractive"):
    """Read all files for one process ("exclusive"/"diffractive") and return
    {(y, pt): [(x_po, dsigma_dy_dpt_dlnxpo), ...]}."""
    name = f"D0_{process}_xpom_{FRAG}_{CHANNEL}_{NUCLEUS}_y*_pt*.dat"
    xpom_dir = f"../output/{CHANNEL}/xpom"
    results = {}
    for filename in sorted(glob.glob(f"{xpom_dir}/{FRAG}/{name}") or glob.glob(f"{xpom_dir}/{name}")):
        y, pt = read_header(filename)
        b_d_list, xpo_list, dsigma_list = read_data_file(filename)
        xpo_groups = group_by_xpo(b_d_list, xpo_list, dsigma_list)
        # The diffractive column is dsigma/dln(x_po); the exclusive one is dsigma/dx_po and gets the factor x_po.
        to_dlnxpo = (lambda xpo: 1.0) if process == "diffractive" else (lambda xpo: xpo)

        points = []
        for xpo, pairs in sorted(xpo_groups.items()):
            b_d_integral = integrate_over_b_d(pairs)
            dsigma = (2*math.pi) * b_d_integral * prefactor(process, pt) \
                     * (2*math.pi) * pt * GEVSQR_TO_MB
            points.append((xpo, to_dlnxpo(xpo) * dsigma))
        results[(y, pt)] = sorted(points)
    return results


def main():
    results = load_results("diffractive")
    if not results:
        sys.exit(f"No files found for NUCLEUS={NUCLEUS}, FRAG_TYPE={FRAG}, CHANNEL={CHANNEL} "
                  f"in ../output/{CHANNEL}/xpom/ -- run ../run_scripts/run_xpom.sh first.")

    PT_TO_PLOT = [
        2.0,
    ]
    filtered_results = {}
    for (y, pt), points in results.items():
        if pt in PT_TO_PLOT:
            filtered_results[(y, pt)] = points
    results = filtered_results

    y_values = sorted(set(y for y, pt in results))
    pt_values = sorted(set(pt for y, pt in results))

    palette = ["#2166ac", "#4393c3", "#92c5de", "#f4a582", "#d6604d", "#b2182b"]
    if len(y_values) <= 4:
        palette = ["#2166ac", "#67a9cf", "#ef8a62", "#b2182b"]
    colors = {y: palette[round(i * (len(palette) - 1) / max(len(y_values) - 1, 1))]
              for i, y in enumerate(y_values)}
    linestyle_cycle = ["-", (0, (8, 2)), "--", "-.", (0, (5, 1.5, 1, 1.5, 1, 1.5)), ":"]
    y_linestyles = {y: linestyle_cycle[i % len(linestyle_cycle)] for i, y in enumerate(y_values)}

    plt.figure(figsize=(8, 7))
    for (y, pt), points in sorted(results.items()):
        xpo_values = [p[0] for p in points]
        dsigma_values = [p[1] for p in points]
        ln_xpo = np.log(xpo_values)
        plt.plot(ln_xpo, dsigma_values, color=colors[y], linestyle=y_linestyles[y], lw=2.0)

    plt.yscale("log")
    plt.ylim(1e-9, 1e5)
    plt.yticks(10.0**np.arange(-9, 6, 2))
    plt.xlim(-10, -2)
    plt.xlabel(r"$\ln(x_{\mathbb{P}})$", labelpad=14)
    plt.ylabel(r"$d\sigma/dy\,dp_{D\perp}\,d\ln(x_{\mathbb{P}})$ [mb/GeV]", labelpad=11)

    ax = plt.gca()
    ax.tick_params(axis='x', which='both', top=False)
    secax = ax.secondary_xaxis('top', functions=(np.exp, np.log))
    secax.set_xscale('log')
    secax.set_xlabel(r"$x_{\mathbb{P}}$", labelpad=23)

    process_label = r"Diffractive$_{\mbox{\fontsize{13}{13}\selectfont TMD}}$"
    pt_label = ", ".join(f"{pt:g}" for pt in pt_values)
    ax.text(0.95, 0.05,
            f"{NUCLEUS}-{NUCLEUS} 5.36 TeV\n{process_label}\n{CHANNEL}, {FRAG}\n"
            f"$p_{{D\\perp}} = {pt_label}$ GeV",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=21, linespacing=1.8)

    y_handles = [Line2D([0], [0], color=colors[y], linestyle=y_linestyles[y], lw=3, label=f"$y={y:g}$")
                 for y in y_values]
    plt.legend(handles=y_handles, loc="upper left", ncol=2, fontsize=21, frameon=False, columnspacing=1.2)

    plt.tight_layout()
    outname = f"../plots/xpom_{CHANNEL}_{NUCLEUS}{NUCLEUS}.pdf"
    plt.savefig(outname, dpi=150)
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
