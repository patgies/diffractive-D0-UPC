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
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

# R_pA at fixed photon energy q+ = 2p+ (z = 0.5), for the charm quark: no photon flux, no fragmentation.
# Pb: run_charm_fixed_qp.sh with DIPOLE_DIR=data/Pb/mve (integrated over b_d); proton: run_charm_fixed_qp.sh.
from D0 import prefactor, nuclear_edge_factor, alphas_run, mc, SIGMA0_MB, GEVSQR_TO_MB, NUCLEUS, NUCLEAR_EDGE
from RpA import A_PB, Y_TO_PLOT

CHARM_DIR = os.environ.get("CHARM_DIR", "../output/charm")
MU_R_FACTORS = [0.5, 1.0, 2.0]   # renormalization-scale variation, mu_R = factor * m_t
# BK posterior samples: member_NNNN/ folders with the layout of output/charm/, each from the Pb dipoles
# data/Pb/bk_posterior/member_NNNN/ and the proton dipole data/proton/bk_posterior/member_NNNN.dat (see README)
BK_DIR = os.environ.get("BK_DIR", "../output/charm_bk")
BK_SAMPLES = "../input/BK/posteriorsamples_100_LOmvefit.dat"
PROCESS_PARTS = {"sum": ["exclusive", "diffractive"], "exclusive": ["exclusive"], "diffractive": ["diffractive"]}


def read_file(filename):
    """Returns (y, pt, exclusive, diffractive), the last three as arrays, without prefactors."""
    with open(filename) as f:
        y = next(float(line.split("=")[-1]) for line in f if "y (fixed)" in line)
    data = np.loadtxt(filename, comments="#")
    return y, data[:, 0], data[:, 1], data[:, 2]


def proton_results(charm_dir=None, sigma0_mb=SIGMA0_MB):
    """{y: {process: (pt, dsigma/d2p in mb/GeV^2)}} for the proton. sigma0_mb is sigma0/2 of the dipole fit."""
    charm_dir = CHARM_DIR if charm_dir is None else charm_dir
    out = {}
    for filename in glob.glob(f"{charm_dir}/charm_fixed_qp_y*.dat"):
        y, pt, excl, diff = read_file(filename)
        out[y] = {proc: (pt, vals * sigma0_mb * np.array([prefactor(proc, p) for p in pt]))
                  for proc, vals in [("exclusive", excl), ("diffractive", diff)]}
    return out


def nucleus_results(charm_dir=None):
    """{y: {process: (pt, dsigma/d2p in mb/GeV^2)}} for the nucleus, integrated over b_d."""
    charm_dir = CHARM_DIR if charm_dir is None else charm_dir
    per_y = {}
    for filename in glob.glob(f"{charm_dir}/{NUCLEUS}/b*/charm_fixed_qp_y*.dat"):
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


def r_pa(aa, pa, y, procs):
    """R_pA at rapidity y for the sum of the processes in procs."""
    return sum(aa[y][p][1] for p in procs) / (A_PB * sum(pa[y][p][1] for p in procs))


def bk_relative_spread():
    """Relative spread (std/mean) of R_pA over the BK posterior samples: for each sample, Pb and the
    proton use the dipoles of the same sample, and the proton its own sigma0/2. The samples come from
    another fit than the central dipole, so only the relative spread is used.
    Returns {process: {y: array}}, or {} if the runs are missing (see BK_DIR)."""
    member_dirs = sorted(glob.glob(f"{BK_DIR}/member_*"))
    if len(member_dirs) < 2 or not os.path.isfile(BK_SAMPLES):
        return {}
    sigma0_samples = np.loadtxt(BK_SAMPLES)[:, 3]   # sigma0/2 in mb, one line per sample
    ratios = {}
    for member_dir in member_dirs:
        index = int(member_dir[-4:])
        aa, pa = nucleus_results(member_dir), proton_results(member_dir, sigma0_samples[index])
        add_edge(aa, pa, member_dir, sigma0_samples[index])
        for process, procs in PROCESS_PARTS.items():
            for y in set(aa) & set(pa):
                ratios.setdefault(process, {}).setdefault(y, []).append(r_pa(aa, pa, y, procs))
    return {process: {y: np.std(values, axis=0) / np.mean(values, axis=0) for y, values in by_y.items()}
            for process, by_y in ratios.items()}


def add_edge(aa, pa, charm_dir, sigma0_mb):
    """Adds the dilute edge of the nucleus to the nucleus results aa, in place: beyond the last b_d of
    the dipole files the cross section is nuclear_edge_factor (D0.py) times the proton one."""
    if not NUCLEAR_EDGE:
        return aa
    b_max = max(float(re.search(r"/b([0-9.]+)/", f).group(1))
                for f in glob.glob(f"{charm_dir}/{NUCLEUS}/b*/charm_fixed_qp_y*.dat"))
    factor = nuclear_edge_factor(b_max, sigma0_mb)
    for y in set(aa) & set(pa):
        for proc in aa[y]:
            pt, values = aa[y][proc]
            aa[y][proc] = (pt, values + factor * np.interp(pt, *pa[y][proc]))
    return aa


def make_plot(process, aa, pa, bk_spread):
    y_values = sorted(y for y in set(aa) & set(pa) if y in Y_TO_PLOT)
    if not y_values:
        sys.exit(f"No common rapidities for Pb and proton in {CHARM_DIR} (run ../run_scripts/run_charm_fixed_qp.sh).")

    y_colors = {0.0: "#2166ac", 0.5: "#67a9cf", 1.0: "#ef8a62", 1.5: "#d6604d", 2.0: "#b2182b"}
    linestyle_cycle = ["-", (0, (8, 2)), "--", "-.", (0, (5, 1.5, 1, 1.5, 1, 1.5)), ":"]
    linestyles = {y: linestyle_cycle[i % len(linestyle_cycle)] for i, y in enumerate(y_values)}

    plt.figure(figsize=(7.5, 6.5))
    for y in y_values:
        procs = PROCESS_PARTS[process]
        pt = aa[y][procs[0]][0]
        central = r_pa(aa, pa, y, procs)
        plt.plot(pt, central, color=y_colors.get(y, "0.3"), linestyle=linestyles[y], linewidth=2)

        # band: BK samples (relative spread), plus the mu_R variation for the sum, in quadrature
        half_width = np.zeros_like(central)
        if y in bk_spread.get(process, {}):
            half_width = bk_spread[process][y] * central
        if process == "sum":
            # mu_R band: alpha_s only multiplies the diffractive part, the same in Pb and in the proton,
            # so it changes the weight of the two processes (no band for each process alone)
            ratios = []
            for mu_r in MU_R_FACTORS:
                rescale = np.array([alphas_run(mu_r * math.sqrt(p**2 + mc**2)) / alphas_run(math.sqrt(p**2 + mc**2))
                                    for p in pt])
                ratios.append((aa[y]["exclusive"][1] + rescale * aa[y]["diffractive"][1])
                              / (A_PB * (pa[y]["exclusive"][1] + rescale * pa[y]["diffractive"][1])))
            mu_r_half_width = (np.maximum.reduce(ratios) - np.minimum.reduce(ratios)) / 2
            half_width = np.sqrt(half_width**2 + mu_r_half_width**2)
        if half_width.any():
            plt.fill_between(pt, central - half_width, central + half_width,
                             color=y_colors.get(y, "0.3"), alpha=0.25, linewidth=0)

    ax = plt.gca()
    plt.xlim(0, 12)
    plt.xticks([0, 2, 4, 6, 8, 10, 12])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    plt.xlabel(r"$p_\perp$ [GeV]", labelpad=15)
    plt.ylabel(r"$R_{pA}$", labelpad=10)

    label = {"sum": "Diffractive", "exclusive": "Exclusive",
             "diffractive": r"Diffractive$_{\mbox{\fontsize{14}{14}\selectfont TMD}}$"}[process]
    # (label top, legend anchor): placed where the curves of each process leave room
    (label_x, label_top, label_ha), legend_loc, legend_anchor, ncol = {
        "sum": ((0.08, 0.9, "left"), "lower right", (0.97, 0.03), 1),
        "exclusive": ((0.08, 0.9, "left"), "lower right", (0.97, 0.03), 1),
        "diffractive": ((0.95, 0.78, "right"), "lower right", (0.97, 0.01), 1)}[process]
    if process == "sum":
        plt.ylim(0.4, 1.1)
    ax.text(label_x, label_top, label + "\n" + r"$q^+ = 2p^+$", transform=ax.transAxes, ha=label_ha, va="top",
            fontsize=23, linespacing=1.6)
    handles = [Line2D([0], [0], color=y_colors.get(y, "0.3"), linestyle=linestyles[y], linewidth=2, label=f"$y={y:g}$")
               for y in y_values]
    band_parts = (["Dipole fit"] if bk_spread.get(process) else []) + ([r"$\mu_R$"] if process == "sum" else [])
    if band_parts:
        handles.append(Patch(facecolor="0.3", alpha=0.25, edgecolor="none", label=r" $\oplus$ ".join(band_parts)))
    plt.legend(handles=handles, loc=legend_loc, bbox_to_anchor=legend_anchor, ncol=ncol, fontsize=20, frameon=False)

    plt.tight_layout()
    suffix = "" if process == "sum" else f"_{process}"
    outname = f"../plots/RpA_fixed_qp{suffix}.pdf"
    plt.savefig(outname, bbox_inches="tight")
    print(f"Saved: {outname}")


def main():
    aa, pa = nucleus_results(), proton_results()
    add_edge(aa, pa, CHARM_DIR, SIGMA0_MB)
    bk_spread = bk_relative_spread()
    for process in ["sum", "exclusive", "diffractive"]:
        make_plot(process, aa, pa, bk_spread)


if __name__ == "__main__":
    main()
