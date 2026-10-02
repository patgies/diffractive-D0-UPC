import math
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import ticker

# Reuses D0.py's data-loading (and picks up its "fancy" rcParams
# styling too, since that's set at import time). Plots the *total* D0 cross
# section (diffractive + exclusive added together) for a few rapidities.
from D0 import (
    load_results, load_hymnd_sum_band, NUCLEUS, CHANNEL, FRAG_TYPES, LINEWIDTH,
)

Y_TO_PLOT = [-1.0, 0.0, 1.0, 2.0]   # the EFF rerun's Pb grid stops at y=2


def main():
    frag = FRAG_TYPES[0]
    results_diff = load_results("diffractive", frag)
    results_excl = load_results("exclusive", frag)

    missing = [y for y in Y_TO_PLOT if y not in results_diff and y not in results_excl]
    if missing:
        sys.exit(f"No data found for rapidities {missing} -- run ../local_workflows/run_many_nucleus.sh first.")

    # Combined "HymnD" theory-uncertainty band for the sum (HymnD-replica
    # and scale-variation uncertainties added in quadrature -- see
    # load_hymnd_sum_band). {} until both run_HymnD_members_roihu.sbatch /
    # run_HymnD_scale_variation.sh have been run and pulled back.
    hymnd_band = load_hymnd_sum_band(frag)

    # Same blue -> red ColorBrewer RdBu steps as D0.py, one per
    # rapidity in increasing y.
    palette = ["#2166ac", "#67a9cf", "#ef8a62", "#b2182b"]
    colors = {y: palette[i % len(palette)] for i, y in enumerate(Y_TO_PLOT)}
    # Linestyle also changes with y, so curves stay distinguishable where
    # they overlap (e.g. y=-1 and y=0) and in black-and-white print.
    linestyle_cycle = ["-", "--", ":", "-."]
    linestyles = {y: linestyle_cycle[i % len(linestyle_cycle)] for i, y in enumerate(Y_TO_PLOT)}

    plt.figure(figsize=(8, 7))
    for y in Y_TO_PLOT:
        band = hymnd_band.get(y)
        if band:
            pt_values, central_values, sigma_total = band
            lower = [v - s for v, s in zip(central_values, sigma_total)]
            upper = [v + s for v, s in zip(central_values, sigma_total)]
            plt.fill_between(pt_values, lower, upper, color=colors[y], alpha=0.25, linewidth=0)
            plt.plot(pt_values, central_values, color=colors[y], linestyle=linestyles[y], lw=LINEWIDTH)
            continue

        diff_points = dict(results_diff.get(y, []))
        excl_points = dict(results_excl.get(y, []))
        pt_values = sorted(set(diff_points) | set(excl_points))
        totals = [diff_points.get(pt, 0.0) + excl_points.get(pt, 0.0) for pt in pt_values]
        plt.plot(pt_values, totals, color=colors[y], linestyle=linestyles[y], lw=LINEWIDTH)

    plt.yscale("log")
    # Fixed range rather than autoscale -- see D0.py's comment:
    # the HymnD band's lower edge crashes toward the 1e-30 floor at a few
    # very-low-pT points (charm-mass threshold effect in the scale
    # variation), which would otherwise stretch the whole axis.
    plt.ylim(1e-8, 1e1)
    plt.xlim(0, 12)   # the pD0 grid of run_many_nucleus.sh ends at 12 GeV
    plt.xlabel(r"$p_{D^0\perp}$ [GeV]", labelpad=14)
    plt.ylabel(r"$d\sigma/dy\,dp_{D^0\perp}$ [mb/GeV]", labelpad=16)
    # in-plot label instead of a title, as in D0.py
    plt.text(0.95, 0.95, f"{NUCLEUS}-{NUCLEUS} 5.36 TeV\n{CHANNEL}, HymnD\nExclusive + diffractive",
             transform=plt.gca().transAxes, ha="right", va="top", fontsize=20, linespacing=1.4)

    # Same numticks fix as D0.py: with a wide log-scale range,
    # LogLocator's default budget silently drops minor ticks entirely, and
    # thins major tick MARKS (not just labels) at every other decade.
    plt.gca().yaxis.set_minor_locator(
        ticker.LogLocator(base=10.0, subs=[2, 3, 4, 5, 6, 7, 8, 9], numticks=100))
    plt.gca().yaxis.set_minor_formatter(ticker.NullFormatter())
    plt.gca().yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=100))
    plt.gca().yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda val, pos: f"$10^{{{round(math.log10(val))}}}$"
                              if round(math.log10(val)) % 2 != 0 else ""))

    # Discrete per-y legend instead of a colorbar (see D0.py).
    y_handles = [Line2D([0], [0], color=colors[y], linestyle=linestyles[y], linewidth=LINEWIDTH, label=f"$y={y:g}$")
                 for y in Y_TO_PLOT]
    plt.legend(handles=y_handles, loc="lower left", bbox_to_anchor=(0.02, 0.0), fontsize=20,
               frameon=False)

    plt.tight_layout()
    outname = f"../plots/D0_sum_{CHANNEL}_{NUCLEUS}{NUCLEUS}.pdf"
    plt.savefig(outname, bbox_inches="tight")
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
