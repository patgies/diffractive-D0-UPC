import glob
import math
import os
import re
import sys
import numpy as np
from scipy.integrate import simpson
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

# R_pA at fixed photon energy q+ = 2p+ (z = 0.5), for the charm quark: no photon flux, no fragmentation.
# Pb: run_charm_fixed_qp.sh with DIPOLE_DIR=data/Pb/mve (integrated over b_d); proton: run_charm_fixed_qp.sh.
from D0 import alphae, mc, e_c, Nc, GEVSQR_TO_MB, NUCLEUS
from RpA import sigma0, A_PB, Y_TO_PLOT
from alphas_running import alphas_run

CHARM_DIR = os.environ.get("CHARM_DIR", "../output/charm")


def prefactor(process, pt):
    """Prefactor without sigma0 (the proton gets sigma0, the nucleus the b_d integral)."""
    if process == "exclusive":
        return alphae * Nc * e_c**2 / (2 * math.pi**2)
    alphas = alphas_run(math.sqrt(pt**2 + mc**2))
    return alphas * alphae * e_c**2 * (Nc**2 - 1) / (8 * math.pi**4)


def read_file(filename):
    """Returns (y, pt, exclusive, diffractive), the last three as arrays, without prefactors."""
    with open(filename) as f:
        y = next(float(line.split("=")[-1]) for line in f if "y (fixed)" in line)
    data = np.loadtxt(filename, comments="#")
    return y, data[:, 0], data[:, 1], data[:, 2]


def proton_results():
    """{y: {process: (pt, dsigma/d2p in mb/GeV^2)}} for the proton."""
    out = {}
    for filename in glob.glob(f"{CHARM_DIR}/charm_fixed_qp_y*.dat"):
        y, pt, excl, diff = read_file(filename)
        out[y] = {proc: (pt, vals * sigma0 * np.array([prefactor(proc, p) for p in pt]))
                  for proc, vals in [("exclusive", excl), ("diffractive", diff)]}
    return out


def nucleus_results():
    """{y: {process: (pt, dsigma/d2p in mb/GeV^2)}} for the nucleus, integrated over b_d."""
    per_y = {}
    for filename in glob.glob(f"{CHARM_DIR}/{NUCLEUS}/b*/charm_fixed_qp_y*.dat"):
        b_d = float(re.search(r"/b([0-9.]+)/", filename).group(1))
        y, pt, excl, diff = read_file(filename)
        per_y.setdefault(y, []).append((b_d, pt, excl, diff))
    out = {}
    for y, samples in per_y.items():
        samples.sort(key=lambda s: s[0])
        if len(samples) < 3:
            continue
        b_d = np.array([s[0] for s in samples])
        pt = samples[0][1]
        out[y] = {}
        for proc, col in [("exclusive", 2), ("diffractive", 3)]:
            values = np.array([s[col] for s in samples])
            b_integral = simpson(b_d[:, None] * values, x=b_d, axis=0)
            out[y][proc] = (pt, 2 * math.pi * b_integral * np.array([prefactor(proc, p) for p in pt]) * GEVSQR_TO_MB)
    return out


def make_plot(process, aa, pa):
    y_values = sorted(y for y in set(aa) & set(pa) if y in Y_TO_PLOT)
    if not y_values:
        sys.exit(f"No common rapidities for Pb and proton in {CHARM_DIR} (run ../run_scripts/run_charm_fixed_qp.sh).")

    y_colors = {0.0: "#2166ac", 0.5: "#67a9cf", 1.0: "#ef8a62", 1.5: "#d6604d", 2.0: "#b2182b"}
    linestyle_cycle = ["-", (0, (8, 2)), "--", "-.", (0, (5, 1.5, 1, 1.5, 1, 1.5)), ":"]
    linestyles = {y: linestyle_cycle[i % len(linestyle_cycle)] for i, y in enumerate(y_values)}

    plt.figure(figsize=(7.5, 6.5))
    for y in y_values:
        procs = ["exclusive", "diffractive"] if process == "sum" else [process]
        pt = aa[y][procs[0]][0]
        num = sum(aa[y][p][1] for p in procs)
        den = sum(pa[y][p][1] for p in procs)
        plt.plot(pt, num / (A_PB * den), color=y_colors.get(y, "0.3"), linestyle=linestyles[y], linewidth=2)

    ax = plt.gca()
    plt.xlim(0, 12)
    plt.xticks([0, 2, 4, 6, 8, 10, 12])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    plt.xlabel(r"$p_\perp$ [GeV]", labelpad=15)
    plt.ylabel(r"$R_{pA}$", labelpad=15)
    plt.title(r"$\gamma + A \to c + A$", fontsize=25, pad=12)

    label = {"sum": "Diffractive", "exclusive": "Exclusive",
             "diffractive": r"Diffractive$_{\mbox{\fontsize{14}{14}\selectfont TMD}}$"}[process]
    # (label top, legend anchor): placed where the curves of each process leave room
    (label_x, label_top, label_ha), legend_loc, legend_anchor, ncol = {
        "sum": ((0.08, 0.9, "left"), "upper right", (0.97, 0.97), 1),   # at the top (y axis extended up)
        "exclusive": ((0.08, 0.9, "left"), "lower right", (0.97, 0.03), 1),
        "diffractive": ((0.95, 0.78, "right"), "lower right", (0.97, 0.01), 1)}[process]
    if process == "sum":
        plt.ylim(0.8, 2.0)
    ax.text(label_x, label_top, label + "\n" + r"$q^+ = 2p^+$", transform=ax.transAxes, ha=label_ha, va="top",
            fontsize=23, linespacing=1.6)
    handles = [Line2D([0], [0], color=y_colors.get(y, "0.3"), linestyle=linestyles[y], linewidth=2, label=f"$y={y:g}$")
               for y in y_values]
    plt.legend(handles=handles, loc=legend_loc, bbox_to_anchor=legend_anchor, ncol=ncol, fontsize=20, frameon=False)

    plt.tight_layout()
    suffix = "" if process == "sum" else f"_{process}"
    outname = f"../plots/RpA_fixed_qp{suffix}.pdf"
    plt.savefig(outname, bbox_inches="tight")
    print(f"Saved: {outname}")


def main():
    aa, pa = nucleus_results(), proton_results()
    for process in ["sum", "exclusive", "diffractive"]:
        make_plot(process, aa, pa)


if __name__ == "__main__":
    main()
