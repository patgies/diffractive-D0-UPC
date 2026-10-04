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
from matplotlib.ticker import MultipleLocator, ScalarFormatter


# Reads the output of run_proton.sh (TARGET=pA).
from D0 import read_rapidity, alphas_run, alphae, mc, e_c, Nc, _scale_combos, SCALE_FACTOR_TAGS
from RpA import sigma0

PROCESS = os.environ.get("PROCESS", "sum")
PPB_DIR = os.environ.get("PPB_DIR", "../output/pPb/central_values")
SCALE_BASE = os.environ.get("SCALE_BASE", "../output/pPb/scale_variation")

PT_BINS = [(0.0, 1.0), (1.0, 2.0), (2.0, 3.0), (3.0, 4.0), (6.0, 7.0), (10.0, 11.0)]
Y_BINS = [(-3.0, -2.0), (-2.0, -1.0), (-1.0, 0.0), (0.0, 1.0), (1.0, 2.0), (2.0, 3.0)]
FRAG_SCHEMES = [
    ("BCFY", "BCFY", "#2166ac", "-"),
    ("KniehlKramer", "KK", "#1b7837", "--"),
    ("HymnD", "HymnD", "#b2182b", ":"),
]
PROCESS_LABELS = {
    "sum": "exclusive + diffractive",
    "exclusive": "exclusive",
    "diffractive": r"diffractive$_{\mbox{\fontsize{13.5}{13.5}\selectfont TMD}}$",
}


class TrimmedFormatter(ScalarFormatter):
    """Tick labels without zeros at the end of a decimal: 0.01 instead of 0.010, 1 instead of 1.0."""
    def _set_format(self):
        super()._set_format()
        self.format = re.sub(r"%1\.\d+f", "%g", self.format)


def proton_prefactor(process, pt, mu_r_factor=1.0):
    """Prefactor for a proton target (sigma0 in mb). Only the diffractive one has alpha_s,
    at the scale mu_r_factor * m_T."""
    if process == "exclusive":
        return alphae * Nc * e_c**2 * sigma0 / (2 * math.pi**2)
    alphas = alphas_run(mu_r_factor * math.sqrt(pt**2 + mc**2))
    return alphas * alphae * e_c**2 * (Nc**2 - 1) * sigma0 / (8 * math.pi**4)


def load_pPb(process, frag, data_dir, mu_r_factor=1.0):
    """Read data_dir/D0_pA_<process>_<frag>*_y*.dat and return {y: [(pt, dsigma_dy_dpt), ...]}."""
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
    """{y: [(pt, cross_section), ...]} for PROCESS and one fragmentation function."""
    if PROCESS != "sum":
        return load_pPb(PROCESS, frag, data_dir, mu_r_factor)
    excl = load_pPb("exclusive", frag, data_dir)
    diff = load_pPb("diffractive", frag, data_dir, mu_r_factor)
    return {y: [(pt, value + dict(diff[y])[pt]) for pt, value in excl[y] if pt in dict(diff[y])]
            for y in excl if y in diff}


def pt_integrate(points, pt_lo, pt_hi):
    pt_values, cs_values = zip(*sorted([(0.0, 0.0)] + list(points)))
    inside = [pt for pt in pt_values if pt_lo < pt < pt_hi]
    x = [pt_lo] + inside + [pt_hi]
    return simpson(np.interp(x, pt_values, cs_values), x=x)


def bin_averages(results):
    """Average of dsigma/(dy dpT) over each (pT, y) bin. Returns {(pt_lo, y_lo): value}."""
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
    """Minimum and maximum of the bin averages over the 7 (mu_F, mu_R) points. {} if there is no data."""
    dirs = {mu_f: f"{PPB_DIR}/{frag}" if tag is None else f"{SCALE_BASE}/{frag}/factor_{tag}"
            for mu_f, tag in SCALE_FACTOR_TAGS.items()}
    if not all(os.path.isdir(d) for d in dirs.values()):
        return {}
    all_values = [bin_averages(get_results(frag, dirs[mu_f], mu_r)) for mu_f, mu_r in _scale_combos()]
    keys = set(all_values[0]).intersection(*all_values[1:])
    return {key: (min(v[key] for v in all_values), max(v[key] for v in all_values)) for key in keys}


def draw_panels(averages, bands, outname, legend_side="right", ylabel_x=0.022, left=0.015, process_text=None):
    """One panel per pT bin with the bin averages {frag: {(pt_lo, y_lo): value}} as histograms in y
    and the bands {frag: {(pt_lo, y_lo): (low, high)}}. Also used by PbPb_bins.py.
    The legend goes in the first panel: upper left, or on the right under the pT label."""
    fig, axes = plt.subplots(3, 2, figsize=(16, 17), sharex=True)
    edges = [y_lo for y_lo, _ in Y_BINS] + [Y_BINS[-1][1]]

    # Texts at the top of a panel: (bins under the text, fraction of the height that is free below it).
    # The vertical range is the smallest one where the histograms stay below these texts.
    on_the_right = lambda y_lo, y_hi: y_hi > -0.4
    on_the_left = lambda y_lo, y_hi: y_lo < -1
    texts = [(on_the_right, 0.77)]   # the pT label
    first_panel_texts = list(texts)
    if legend_side == "right":
        first_panel_texts.append((lambda y_lo, y_hi: y_hi > 0, 0.43))   # legend under the pT label
        if process_text:
            first_panel_texts.append((on_the_left, 0.84))
    else:
        first_panel_texts.append((on_the_left, 0.52))                  # legend in the upper left
        if process_text:
            first_panel_texts.append((on_the_right, 0.72))

    for ax, (pt_lo, pt_hi) in zip(axes.flat, PT_BINS):
        bin_top = {y_bin: 0.0 for y_bin in Y_BINS}   # highest point of each y bin
        for frag, _, color, linestyle in FRAG_SCHEMES:
            values = [averages[frag].get((pt_lo, y_lo), np.nan) for y_lo, _ in Y_BINS]
            if np.all(np.isnan(values)):
                continue
            ax.stairs(values, edges, baseline=None, color=color, linestyle=linestyle, lw=3)
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
        for under_text, free_below in (first_panel_texts if ax is axes.flat[0] else texts):
            y_max = max([y_max] + [top / free_below for y_bin, top in bin_top.items() if under_text(*y_bin)])
        ax.set_ylim(0, y_max)
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.tick_params(labelsize=30, pad=10)
        ax.tick_params(which="both", color="black")
        formatter = TrimmedFormatter(useMathText=True)
        if pt_hi <= 2.0:
            formatter.set_scientific(False)   # decimals in the first two panels
        else:
            formatter.set_powerlimits((-2, 2))
            ax.yaxis.get_offset_text().set_fontsize(30)
        ax.yaxis.set_major_formatter(formatter)
        ax.text(0.95, 0.93, rf"${pt_lo:g} < p_{{D^0\perp}} < {pt_hi:g}$ GeV", transform=ax.transAxes,
                ha="right", va="top", fontsize=25)

    frag_handles = [Line2D([0], [0], color=color, linestyle=linestyle, lw=3, label=label)
                    for _, label, color, linestyle in FRAG_SCHEMES]
    if process_text:   # in the first panel, on the side that the legend does not use
        x, y, ha = (0.05, 0.93, "left") if legend_side == "right" else (0.95, 0.82, "right")
        axes.flat[0].text(x, y, process_text, transform=axes.flat[0].transAxes, ha=ha, va="top", fontsize=25)
    if legend_side == "right":
        axes.flat[0].legend(handles=frag_handles, loc="upper right", bbox_to_anchor=(1.0, 0.87), fontsize=22,
                            frameon=False)
    else:
        axes.flat[0].legend(handles=frag_handles, loc="upper left", fontsize=22, frameon=False)


    for ax in axes[-1]:
        ax.set_xlabel(r"$y$", labelpad=12, fontsize=34)
    fig.supylabel(r"$d\sigma/dy\,dp_{D^0\perp}$ [mb/GeV]", fontsize=34, x=ylabel_x)

    plt.tight_layout(h_pad=1.5, w_pad=3.0, rect=(left, 0, 1, 1))
    plt.savefig(outname, bbox_inches="tight")
    print(f"Saved: {outname}")


def main():
    averages = {frag: bin_averages(get_results(frag, f"{PPB_DIR}/{frag}")) for frag, _, _, _ in FRAG_SCHEMES}
    bands = {frag: scale_band(frag) for frag, _, _, _ in FRAG_SCHEMES}
    if not any(averages.values()):
        sys.exit(f"No data found in {PPB_DIR} -- run TARGET=pA ../run_scripts/run_proton.sh first.")

    suffix = "" if PROCESS == "sum" else f"_{PROCESS}"
    draw_panels(averages, bands, f"../plots/D0_bins_pPb{suffix}.pdf", legend_side="left")


if __name__ == "__main__":
    main()
