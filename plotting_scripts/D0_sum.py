import math
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import ticker

# Total D0 cross section (diffractive + exclusive) for some rapidities,

from D0 import (
    load_results, load_hymnd_sum_band, NUCLEUS, CHANNEL, FRAG_TYPES, LINEWIDTH,
)

Y_TO_PLOT = [-1.0, 0.0, 1.0, 2.0]   


def main():
    frag = FRAG_TYPES[0]
    results_diff = load_results("diffractive", frag)
    results_excl = load_results("exclusive", frag)

    missing = [y for y in Y_TO_PLOT if y not in results_diff and y not in results_excl]
    if missing:
        sys.exit(f"No data found for rapidities {missing} -- run ../run_scripts/run_nucleus.sh first.")

    hymnd_band = load_hymnd_sum_band(frag)

    palette = ["#2166ac", "#67a9cf", "#ef8a62", "#b2182b"]
    colors = {y: palette[i % len(palette)] for i, y in enumerate(Y_TO_PLOT)}
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
    plt.ylim(3e-8, 1.5 * plt.gca().dataLim.y1)   # just above the highest curve
    plt.xlim(0, 12)   
    plt.xlabel(r"$p_{D^0\perp}$ [GeV]", labelpad=14)
    plt.ylabel(r"$d\sigma/dy\,dp_{D^0\perp}$ [mb/GeV]", labelpad=16)
    plt.text(0.95, 0.92, f"{NUCLEUS}-{NUCLEUS} 5.36 TeV\n{CHANNEL}, HymnD",
             transform=plt.gca().transAxes, ha="right", va="top", fontsize=24, linespacing=1.8)

    plt.gca().yaxis.set_minor_locator(
        ticker.LogLocator(base=10.0, subs=[2, 3, 4, 5, 6, 7, 8, 9], numticks=100))
    plt.gca().yaxis.set_minor_formatter(ticker.NullFormatter())
    plt.gca().yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=100))
    plt.gca().yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda val, pos: f"$10^{{{round(math.log10(val))}}}$"
                              if round(math.log10(val)) % 2 != 0 else ""))

    y_handles = [Line2D([0], [0], color=colors[y], linestyle=linestyles[y], linewidth=3, label=f"$y={y:g}$")
                 for y in Y_TO_PLOT]
    plt.legend(handles=y_handles, loc="lower left", bbox_to_anchor=(0.02, 0.0), fontsize=22,
               frameon=False)

    plt.tight_layout()
    outname = f"../plots/D0_sum_{CHANNEL}_{NUCLEUS}{NUCLEUS}.pdf"
    plt.savefig(outname, bbox_inches="tight")
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
