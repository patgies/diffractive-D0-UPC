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
from matplotlib.ticker import MultipleLocator, ScalarFormatter


# Input: run_proton.sh (TARGET=pA).
from D0 import read_rapidity, proton_prefactor, sum_with_limit, _scale_combos, SCALE_FACTOR_TAGS

PROCESS = os.environ.get("PROCESS", "sum")
PPB_DIR = os.environ.get("PPB_DIR", "../output/pPb/central_values")
SCALE_BASE = os.environ.get("SCALE_BASE", "../output/pPb/scale_variation")

PT_BINS = [(1.0, 2.0), (2.0, 3.0), (3.0, 4.0), (6.0, 7.0), (10.0, 11.0)]
Y_BINS = [(-3.0, -2.0), (-2.0, -1.0), (-1.0, 0.0), (0.0, 1.0), (1.0, 2.0), (2.0, 3.0)]
FRAG_SCHEMES = [
    ("BCFY", "BCFY", "#2166ac", "-"),
    ("KniehlKramer", "KK", "#1b7837", "--"),
    ("HymnD", "HymnD", "#b2182b", ":"),
]
PANEL_PROCESS_LABELS = {
    "sum": "Diffractive",
    "exclusive": "Exclusive",
    "diffractive": r"Diffractive$_{\mbox{\fontsize{13.5}{13.5}\selectfont TMD}}$",
}


class TrimmedFormatter(ScalarFormatter):
    """Tick labels without trailing zeros."""
    def _set_format(self):
        super()._set_format()
        self.format = re.sub(r"%1\.\d+f", "%g", self.format)


def load_pPb(process, frag, data_dir, mu_r_factor=1.0):
    """{y: [(pt, dsigma)]} from D0_pA files."""
    results = {}
    for filename in sorted(glob.glob(f"{data_dir}/D0_pA_{process}_{frag}*_y*.dat")):
        y = read_rapidity(filename)
        points = {}
        with open(filename) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                pt, raw = line.split()
                pt, raw = float(pt), float(raw)
                points[pt] = raw * proton_prefactor(process, pt, mu_r_factor) * (2 * math.pi) * pt
        results[y] = sorted(points.items())
    return results


def get_results(frag, data_dir, mu_r_factor=1.0):
    """{y: [(pt, cross_section), ...]} for one fragmentation function."""
    if PROCESS != "sum":
        return load_pPb(PROCESS, frag, data_dir, mu_r_factor)
    excl = load_pPb("exclusive", frag, data_dir)
    diff = load_pPb("diffractive", frag, data_dir, mu_r_factor)
    return {y: sum_with_limit(dict(excl[y]), dict(diff[y])) for y in excl if y in diff}


def pt_integrate(points, pt_lo, pt_hi):
    pt_values, cs_values = zip(*sorted([(0.0, 0.0)] + list(points)))
    inside = [pt for pt in pt_values if pt_lo < pt < pt_hi]
    x = [pt_lo] + inside + [pt_hi]
    return simpson(np.interp(x, pt_values, cs_values), x=x)


def bin_averages(results):
    """Average over each (pT, y) bin. Returns {(pt_lo, y_lo): value}."""
    out = {}
    for pt_lo, pt_hi in PT_BINS:
        per_y = {y: pt_integrate(points, pt_lo, pt_hi) for y, points in results.items()}
        for y_lo, y_hi in Y_BINS:
            y_list = sorted(y for y in per_y if y_lo - 1e-9 <= y <= y_hi + 1e-9)
            if len(y_list) < 3:
                continue
            integral = simpson([per_y[y] for y in y_list], x=y_list)
            out[(pt_lo, y_lo)] = integral / (y_hi - y_lo) / (pt_hi - pt_lo)
    return out


def scale_band(frag):
    """Min and max of the bin averages over the 7 scale points. {} if no data."""
    dirs = {mu_f: f"{PPB_DIR}/{frag}" if tag is None else f"{SCALE_BASE}/{frag}/factor_{tag}"
            for mu_f, tag in SCALE_FACTOR_TAGS.items()}
    if not all(os.path.isdir(d) for d in dirs.values()):
        return {}
    all_values = [bin_averages(get_results(frag, dirs[mu_f], mu_r)) for mu_f, mu_r in _scale_combos()]
    keys = set(all_values[0]).intersection(*all_values[1:])
    return {key: (min(v[key] for v in all_values), max(v[key] for v in all_values)) for key in keys}


def draw_panels(averages, bands, outname, ylabel_x=0.034, left=0.015, flux_text=None, second_panel_text=None):
    """One panel per pT bin: histograms in y and scale bands. Also used by PbPb_bins.py."""
    fig, axes = plt.subplots(3, 2, figsize=(16, 17), sharex=True)

    # Texts at the top of a panel: (bins under it, free fraction below it).
    on_the_right = lambda y_lo, y_hi: y_hi > -0.4
    on_the_left = lambda y_lo, y_hi: y_lo < -1
    has_band = any(bands[frag] for frag in bands)
    texts = [(on_the_right, 0.77)]
    # Legend in the first free panel, or in the upper left of the first one.
    n_free = axes.size - len(PT_BINS)
    legend_panel = n_free > 0
    bin_axes = axes.flat[:len(PT_BINS)]
    first_panel_texts = texts if legend_panel else [(on_the_right, 0.77), (on_the_left, 0.42 if has_band else 0.6),
                                                    (lambda y_lo, y_hi: y_lo < 1, 0.5 if has_band else 1.0)]

    for ax, (pt_lo, pt_hi) in zip(bin_axes, PT_BINS):
        bin_top = {y_bin: 0.0 for y_bin in Y_BINS}
        for frag, _, color, linestyle in FRAG_SCHEMES:
            values = [averages[frag].get((pt_lo, y_lo), np.nan) for y_lo, _ in Y_BINS]
            if np.all(np.isnan(values)):
                continue
            for value, (y_lo, y_hi) in zip(values, Y_BINS):
                if not np.isnan(value):
                    ax.plot([y_lo, y_hi], [value, value], color=color, linestyle=linestyle, lw=3)
            for value, y_bin in zip(values, Y_BINS):
                if not np.isnan(value):
                    bin_top[y_bin] = max(bin_top[y_bin], value)
            for y_lo, y_hi in Y_BINS:
                if (pt_lo, y_lo) in bands[frag]:
                    low, high = bands[frag][(pt_lo, y_lo)]
                    ax.fill_between([y_lo, y_hi], low, high, color=color, alpha=0.25, linewidth=0)
                    bin_top[(y_lo, y_hi)] = max(bin_top[(y_lo, y_hi)], high)

        ax.set_xlim(-3, 3)
        y_max = 1.12 * max(bin_top.values())
        panel_texts = first_panel_texts if ax is bin_axes[0] else texts
        if flux_text and ax is axes.flat[0] and not legend_panel:
            panel_texts = panel_texts + [(on_the_right, 0.66)]
        if second_panel_text and ax is axes.flat[1] and not legend_panel:
            panel_texts = panel_texts + [(on_the_right, 0.66)]
        for under_text, free_below in panel_texts:
            y_max = max([y_max] + [top / free_below for y_bin, top in bin_top.items() if under_text(*y_bin)])
        ax.set_ylim(0, y_max)
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.tick_params(labelsize=30, pad=10)
        formatter = TrimmedFormatter(useMathText=True)
        if pt_hi <= 2.0:
            formatter.set_scientific(False)
        else:
            formatter.set_powerlimits((-2, 2))
            ax.yaxis.get_offset_text().set(fontsize=30, color="black")
            ax.yaxis.OFFSETTEXTPAD = 8
        ax.yaxis.set_major_formatter(formatter)
        ax.text(0.95, 0.93, rf"$p_{{D\perp}} \in ({pt_lo:g}, {pt_hi:g})$ GeV", transform=ax.transAxes,
                ha="right", va="top", fontsize=30)

    frag_handles = [Line2D([0], [0], color=color, linestyle=linestyle, lw=3, label=label)
                    for _, label, color, linestyle in FRAG_SCHEMES]
    if has_band:
        frag_handles.append(Patch(facecolor="0.3", alpha=0.25, edgecolor="none",
                                  label=r"$\mu_F, \mu_R \in [0.5, 2]\, m_T^D$"))
    if flux_text and not legend_panel:
        axes.flat[0].text(0.95, 0.79, flux_text, transform=axes.flat[0].transAxes, ha="right", va="top", fontsize=30)
    if second_panel_text and not legend_panel:
        axes.flat[1].text(0.95, 0.79, second_panel_text, transform=axes.flat[1].transAxes, ha="right", va="top",
                          fontsize=30)
    if not legend_panel:
        axes.flat[0].legend(handles=frag_handles, loc="upper left", fontsize=28, frameon=False)

    # empty panels hidden, the first one holds the legend
    for i, ax in enumerate(axes.flat[len(PT_BINS):], start=len(PT_BINS)):
        ax.axis("off")
        if i == len(PT_BINS):
            text_handles = [Patch(visible=False, label=t) for t in (second_panel_text, flux_text) if t]
            legend = ax.legend(handles=text_handles + frag_handles, loc="center", fontsize=30,
                               frameon=True, fancybox=False, edgecolor="0.3", borderpad=0.8)
            legend.get_frame().set_linewidth(1.2)
        axes.flat[i - axes.shape[1]].tick_params(labelbottom=True)
        axes.flat[i - axes.shape[1]].set_xlabel(r"$y$", labelpad=12, fontsize=34)
    for ax in axes[-1][:axes.shape[1] - n_free]:
        ax.set_xlabel(r"$y$", labelpad=12, fontsize=34)
    fig.supylabel(r"$d\sigma/dy\,dp_{D\perp}$ [mb/GeV]", fontsize=34, x=ylabel_x)

    plt.tight_layout(h_pad=0.3, w_pad=3.0, rect=(left, 0, 1, 1))
    plt.savefig(outname, bbox_inches="tight")
    print(f"Saved: {outname}")


def main():
    averages = {frag: bin_averages(get_results(frag, f"{PPB_DIR}/{frag}")) for frag, _, _, _ in FRAG_SCHEMES}
    bands = {frag: scale_band(frag) for frag, _, _, _ in FRAG_SCHEMES}
    if not any(averages.values()):
        sys.exit(f"No data found in {PPB_DIR} -- run TARGET=pA ../run_scripts/run_proton.sh first.")

    suffix = "" if PROCESS == "sum" else f"_{PROCESS}"
    draw_panels(averages, bands, f"../plots/D0_bins_pPb{suffix}.pdf",
                second_panel_text=f"{PANEL_PROCESS_LABELS[PROCESS]} p-Pb")


if __name__ == "__main__":
    main()
