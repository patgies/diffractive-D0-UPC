import glob
import math
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
from matplotlib.ticker import FuncFormatter

from D0 import (
    read_rapidity, alphae, mc, e_c, Nc, load_results, _summed_results,
    load_bk_band, load_bk_sum_band, CHANNEL, CENTRAL_DIR,
)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))   # alphas_running.py
from alphas_running import alphas_run

sigma0 = 16.36


def proton_prefactor(process, pt):
    """Prefactor for a proton target. sigma0 is already in mb (no GEVSQR_TO_MB)."""
    if process == "exclusive":
        return alphae * Nc * e_c**2 * sigma0 / (2 * math.pi**2)
    elif process == "diffractive":
        alphas = alphas_run(math.sqrt(pt**2 + mc**2))
        return alphas * alphae * e_c**2 * (Nc**2 - 1) * sigma0 / (8 * math.pi**4)
    raise ValueError(f"unknown process {process}")

# R_pA = dsigma_AA / (A * dsigma_pA), both with the same nuclear photon flux. Numerator:
# Pb+Pb sum (D0.py); denominator: run_proton.sh (TARGET=AA) output, no b_d integral.

A_PB = 208
Y_TO_PLOT = [0.0, 0.5, 1.0, 1.5, 2.0]

PROTON_DIR = os.environ.get("PROTON_DIR", f"../output/{CHANNEL}/proton_baseline")


def load_proton_baseline(process, frag="HymnD"):
    """Read PROTON_DIR/<frag>/D0_proton_baseline_<process>_<frag>_<channel>_y*.dat (no b_d column)
    and return {y: [(pt, dsigma_pA_dy_dpt), ...]}."""
    results = {}
    name = f"D0_proton_baseline_{process}_{frag}_{CHANNEL}_y*.dat"
    for filename in sorted(glob.glob(f"{PROTON_DIR}/{frag}/{name}") or glob.glob(f"{PROTON_DIR}/{name}")):
        y = read_rapidity(filename)
        points = {}
        with open(filename) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                pt, raw = line.split()
                pt, raw = float(pt), float(raw)
                points[pt] = raw * proton_prefactor(process, pt) * (2 * math.pi) * pt
        results[y] = sorted(points.items())
    return results


def load_proton_baseline_sum(frag="HymnD"):
    """Same as load_proton_baseline, but summing exclusive+diffractive."""
    excl = load_proton_baseline("exclusive", frag)
    diff = load_proton_baseline("diffractive", frag)
    results = {}
    for y in set(excl) | set(diff):
        excl_points = dict(excl.get(y, []))
        diff_points = dict(diff.get(y, []))
        pt_values = sorted(set(excl_points) | set(diff_points))
        results[y] = [(pt, excl_points.get(pt, 0.0) + diff_points.get(pt, 0.0)) for pt in pt_values]
    return results


def make_plot(process):
    """process is "exclusive"/"diffractive" (-> load_results/load_proton_baseline)
    or "sum" (-> _summed_results/load_proton_baseline_sum)."""
    if process == "sum":
        aa_results = _summed_results(CENTRAL_DIR, "HymnD")
        pa_results = load_proton_baseline_sum("HymnD")
        bk_band = load_bk_sum_band("HymnD")
        label = "exclusive + diffractive"
    else:
        aa_results = load_results(process, "HymnD")
        pa_results = load_proton_baseline(process, "HymnD")
        bk_band = load_bk_band(process, "HymnD")
        label = process

    y_values = sorted(y for y in set(aa_results) & set(pa_results) if y in Y_TO_PLOT)
    if not y_values:
        sys.exit(f"No overlapping rapidities between Pb+Pb and proton-baseline data for {process} -- "
                  f"check {CENTRAL_DIR} and {PROTON_DIR} (run ../run_scripts/run_proton.sh).")

    y_colors = {-2.0: "#053061", -1.0: "#053061", -0.5: "#2166ac", 0.0: "#2166ac", 0.5: "#67a9cf",
                1.0: "#ef8a62", 1.5: "#d6604d", 2.0: "#b2182b"}
    colors = {y: y_colors.get(y, "0.3") for y in y_values}

    linestyle_cycle = ["-", (0, (8, 2)), "--", "-.", (0, (5, 1.5, 1, 1.5, 1, 1.5)), ":"]
    y_linestyles = {y: linestyle_cycle[i % len(linestyle_cycle)] for i, y in enumerate(y_values)}

    plt.figure(figsize=(7.5, 6.5))
    for y in y_values:
        aa_points = dict(aa_results[y])
        pa_points = dict(pa_results[y])
        pt_values = sorted(set(aa_points) & set(pa_points))
        r_pa = [aa_points[pt] / (A_PB * pa_points[pt]) for pt in pt_values]
        plt.plot(pt_values, r_pa, color=colors[y], linestyle=y_linestyles[y], linewidth=2)

        if y in bk_band:
            bk_pt, bk_mean, bk_std = bk_band[y]
            rel_bk = {pt: s / m for pt, m, s in zip(bk_pt, bk_mean, bk_std) if m > 0}
            band_pt = [pt for pt in pt_values if pt in rel_bk]
            band_r = [aa_points[pt] / (A_PB * pa_points[pt]) for pt in band_pt]
            band_rel = [rel_bk[pt] for pt in band_pt]
            lower = [r * (1 - rel) for r, rel in zip(band_r, band_rel)]
            upper = [r * (1 + rel) for r, rel in zip(band_r, band_rel)]
            plt.fill_between(band_pt, lower, upper, color=colors[y], alpha=0.25, linewidth=0)

    plt.axhline(1.0, color="gray", linestyle=":", linewidth=1)
    plt.xlim(0, 12)
    plt.xticks([0, 2, 4, 6, 8, 10, 12])
    if process == "diffractive":
        plt.ylim(0.4, 1.15)
    else:
        plt.ylim(bottom=0.85)
    plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))

    plt.xlabel(r"$p_{D^0\perp}$ [GeV]", labelpad=15)
    plt.ylabel(r"$R_{pA}$", labelpad=15)

    y_handles = [Line2D([0], [0], color=colors[y], linestyle=y_linestyles[y], linewidth=2, label=f"$y={y:g}$")
                 for y in y_values]
    ax = plt.gca()
    if process == "exclusive":
        ax.text(0.08, 0.93, "Exclusive", transform=ax.transAxes, ha="left", va="top", fontsize=23)
    if process == "sum":
        ax.text(0.95, 0.93, "Diffractive", transform=ax.transAxes, ha="right", va="top", fontsize=23)
    if process == "diffractive":
        ax.text(0.95, 0.93, r"Diffractive$_{\mbox{\fontsize{14}{14}\selectfont TMD}}$", transform=ax.transAxes, ha="right", va="top", fontsize=23)
        plt.legend(handles=y_handles, loc="lower right", bbox_to_anchor=(0.97, 0.03), fontsize=20, frameon=False)
    else:
        plt.legend(handles=y_handles, loc="lower right", bbox_to_anchor=(0.97, 1.16 if process == "sum" else 1.04),
                   fontsize=20, frameon=False, bbox_transform=blended_transform_factory(ax.transAxes, ax.transData))

    plt.tight_layout()
    suffix = "" if process == "sum" else f"_{process}"
    outname = f"../plots/RpA_{CHANNEL}{suffix}.pdf"
    plt.savefig(outname, bbox_inches="tight")
    print(f"Saved: {outname}")


def main():
    make_plot("sum")
    make_plot("exclusive")
    make_plot("diffractive")


if __name__ == "__main__":
    main()
