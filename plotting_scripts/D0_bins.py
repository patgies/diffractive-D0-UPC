import os
import sys
import numpy as np
from scipy.integrate import simpson
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

# D0 cross section averaged over the CMS (pT, y) bins, for the three fragmentation functions,
# with the scale-variation band. 
from D0 import (
    load_results, _summed_results, _scale_combos, SCALE_FACTOR_TAGS,
    NUCLEUS, CHANNEL, CENTRAL_DIR,
)

PROCESS = os.environ.get("PROCESS", "sum")
SCALE_BASE = os.environ.get("SCALE_BASE", f"../output/{CHANNEL}/scale_variation")

PT_BINS_Y_BINS = [
    (2.0, 5.0, [(-1.0, 1.0)]),
    (5.0, 8.0, [(-2.0, -1.0), (-1.0, 0.0), (0.0, 1.0), (1.0, 2.0)]),
    (8.0, 12.0, [(-2.0, -1.0), (-1.0, 0.0), (0.0, 1.0), (1.0, 2.0)]),
]
COLORS = {2.0: "#2166ac", 5.0: "#e08a00", 8.0: "#b2182b"}
FRAG_SCHEMES = [
    ("BCFY", "BCFY", "+", -0.09),
    ("KniehlKramer", "Kniehl-Kramer", "o", 0.0),
    ("HymnD", "HymnD", ".", 0.09),
]
PROCESS_LABELS = {
    "sum": "exclusive + diffractive",
    "exclusive": "exclusive",
    "diffractive": r"diffractive$_{\mbox{\fontsize{13.5}{13.5}\selectfont TMD}}$",
}


def get_results(frag, data_dir, mu_r_factor=1.0):
    """{y: [(pt, cross_section), ...]} for PROCESS and one fragmentation function."""
    if PROCESS == "sum":
        return _summed_results(data_dir, frag, mu_r_factor=mu_r_factor)
    return load_results(PROCESS, frag, data_dir=data_dir, mu_r_factor=mu_r_factor)


def pt_integrate(points, pt_lo, pt_hi):
    pt_values, cs_values = zip(*sorted(points))
    inside = [pt for pt in pt_values if pt_lo < pt < pt_hi]
    x = [pt_lo] + inside + [pt_hi]
    return simpson(np.interp(x, pt_values, cs_values), x=x)


def bin_averages(results):
    """Average of dsigma/(dy dpT) over each (pT, y) bin. Returns {(pt_lo, y_lo): value}."""
    out = {}
    for pt_lo, pt_hi, y_bins in PT_BINS_Y_BINS:
        per_y = {y: pt_integrate(points, pt_lo, pt_hi) for y, points in results.items()}
        for y_lo, y_hi in y_bins:
            y_list = sorted(y for y in per_y if y_lo - 1e-9 <= y <= y_hi + 1e-9)
            if len(y_list) < 3:
                continue
            integral = simpson([per_y[y] for y in y_list], x=y_list)
            out[(pt_lo, y_lo)] = integral / (y_hi - y_lo) / (pt_hi - pt_lo)
    return out


def scale_band(frag):
    """Minimum and maximum of the bin averages over the 7 (mu_F, mu_R) points. {} if there is no data."""
    dirs = {mu_f: CENTRAL_DIR if tag is None else f"{SCALE_BASE}/{frag}/factor_{tag}"
            for mu_f, tag in SCALE_FACTOR_TAGS.items()}
    if not all(os.path.isdir(d) for d in dirs.values()):
        return {}
    all_values = [bin_averages(get_results(frag, dirs[mu_f], mu_r)) for mu_f, mu_r in _scale_combos()]
    keys = set(all_values[0]).intersection(*all_values[1:])
    return {key: (min(v[key] for v in all_values), max(v[key] for v in all_values)) for key in keys}


def main():
    plt.figure(figsize=(8, 7))
    ax = plt.gca()
    found = False

    for frag, frag_label, marker, dx in FRAG_SCHEMES:
        central = bin_averages(get_results(frag, CENTRAL_DIR))
        band = scale_band(frag)
        for pt_lo, pt_hi, y_bins in PT_BINS_Y_BINS:
            color = COLORS[pt_lo]
            for y_lo, y_hi in y_bins:
                if (pt_lo, y_lo) not in central:
                    continue
                found = True
                x = 0.5 * (y_lo + y_hi) + dx
                value = central[(pt_lo, y_lo)]
                if (pt_lo, y_lo) in band:
                    low, high = band[(pt_lo, y_lo)]
                    ax.plot([x, x], [low, high], color=color, lw=2, zorder=1)
                    for edge in (low, high):
                        ax.plot([x - 0.035, x + 0.035], [edge, edge], color=color, lw=2, zorder=1)
                if marker == "+":
                    ax.plot(x, value, "+", color=color, markersize=13, markeredgewidth=2.2, zorder=3)
                elif marker == "o":
                    ax.plot(x, value, "o", color=color, markerfacecolor="white", markersize=9,
                            markeredgewidth=2, zorder=3)
                else:
                    ax.plot(x, value, "o", color=color, markersize=8, zorder=3)

    if not found:
        sys.exit(f"No data found in {CENTRAL_DIR}.")

    frag_handles = [
        Line2D([0], [0], color="0.2", marker="+", linestyle="", markersize=13, markeredgewidth=2.2, label="BCFY"),
        Line2D([0], [0], color="0.2", marker="o", linestyle="", markerfacecolor="white", markersize=9,
               markeredgewidth=2, label="Kniehl-Kramer"),
        Line2D([0], [0], color="0.2", marker="o", linestyle="", markersize=8, label="HymnD"),
    ]
    frag_legend = plt.legend(handles=frag_handles, loc="lower left", bbox_to_anchor=(0.46, 0.0),
                             fontsize=20, frameon=False, handletextpad=0.3)
    ax.add_artist(frag_legend)

    pt_handles = [Line2D([0], [0], color=COLORS[pt_lo], linewidth=3,
                         label=rf"${pt_lo:g} < p_{{D^0\perp}} < {pt_hi:g}$ GeV")
                  for pt_lo, pt_hi, _ in PT_BINS_Y_BINS]
    plt.legend(handles=pt_handles, loc="lower left", bbox_to_anchor=(0.0, 0.0), fontsize=20,
               frameon=False, handlelength=1.2)

    plt.yscale("log")
    plt.xlim(-2.4, 2.4)
    plt.ylim(3e-9, 1e0)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    plt.xlabel(r"$y$", labelpad=10)
    plt.ylabel(r"$d\sigma/dy\,dp_{D^0\perp}$ [mb/GeV]", labelpad=16)
    plt.title(f"{NUCLEUS}-{NUCLEUS} 5.36 TeV, {CHANNEL}, {PROCESS_LABELS[PROCESS]}", pad=15, fontsize=24)

    plt.tight_layout()
    suffix = "" if PROCESS == "sum" else f"_{PROCESS}"
    outname = f"../plots/D0_bins_{CHANNEL}_{NUCLEUS}{NUCLEUS}{suffix}.pdf"
    plt.savefig(outname, bbox_inches="tight")
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
