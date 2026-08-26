import math
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import ticker

# Reuses cross_section.py's data-loading (and picks up its "fancy" rcParams
# styling too, since that's set at import time). Plots the *total* D0 cross
# section (diffractive + exclusive added together) for a few rapidities.
from cross_section import (
    load_results, load_hymnd_sum_band, NUCLEUS, CHANNEL, FRAG_TYPES,
)

Y_TO_PLOT = [-1.0, 0.0, 1.0, 2.0, 3.0]


def main():
    frag = FRAG_TYPES[0]
    results_diff = load_results("diffractive", frag)
    results_excl = load_results("exclusive", frag)

    missing = [y for y in Y_TO_PLOT if y not in results_diff and y not in results_excl]
    if missing:
        sys.exit(f"No data found for rapidities {missing} -- run ../run_many_nucleus.sh first.")

    # Combined "HymnD" theory-uncertainty band for the sum (LHAPDF-replica
    # and scale-variation uncertainties added in quadrature -- see
    # load_hymnd_sum_band). {} until both run_lhapdf_members_roihu.sbatch /
    # run_lhapdf_scale_variation.sh have been run and pulled back.
    hymnd_band = load_hymnd_sum_band(frag)

    # Same ColorBrewer "Blues" scale as cross_section.py.
    blue_ramp = [
        "#deebf7", "#c6dbef", "#9ecae1", "#6baed6",
        "#4292c6", "#2171b5", "#08519c", "#08306b",
    ]
    colors = {}
    for i, y in enumerate(Y_TO_PLOT):
        step = round(i * (len(blue_ramp) - 1) / max(len(Y_TO_PLOT) - 1, 1))
        colors[y] = blue_ramp[step]

    # Color AND linestyle both cycle with y (see fixedW_xpom_spectrum.py) --
    # makes adjacent/crossing curves easier to tell apart than color alone.
    linestyle_cycle = ['-', '--', ':', '-.']
    y_linestyles = {}
    for i, y in enumerate(Y_TO_PLOT):
        y_linestyles[y] = linestyle_cycle[i % len(linestyle_cycle)]

    plt.figure(figsize=(7.5, 6.5))
    for y in Y_TO_PLOT:
        band = hymnd_band.get(y)
        if band:
            pt_values, central_values, sigma_total = band
            lower = [v - s for v, s in zip(central_values, sigma_total)]
            upper = [v + s for v, s in zip(central_values, sigma_total)]
            plt.fill_between(pt_values, lower, upper, color=colors[y], alpha=0.25, linewidth=0)
            plt.plot(pt_values, central_values, color=colors[y], linestyle=y_linestyles[y], linewidth=2)
            continue

        diff_points = dict(results_diff.get(y, []))
        excl_points = dict(results_excl.get(y, []))
        pt_values = sorted(set(diff_points) | set(excl_points))
        totals = [diff_points.get(pt, 0.0) + excl_points.get(pt, 0.0) for pt in pt_values]
        plt.plot(pt_values, totals, color=colors[y], linestyle=y_linestyles[y], linewidth=2)

    plt.yscale("log")
    # Fixed range rather than autoscale -- see cross_section.py's comment:
    # the HymnD band's lower edge crashes toward the 1e-30 floor at a few
    # very-low-pT points (charm-mass threshold effect in the scale
    # variation), which would otherwise stretch the whole axis.
    plt.ylim(1e-10, 1e2)
    plt.xlabel(r"$p_{D^0\perp}$ [GeV]", labelpad=15)
    plt.ylabel(r"$d\sigma/dy\,dp_{D^0\perp}$ [mb/GeV]", labelpad=15)
    plt.title(f"$D^0$ photoproduction (diff. + excl.), {NUCLEUS}+{NUCLEUS} UPC ({CHANNEL}, HymnD)", pad=15)

    # Same numticks fix as cross_section.py: with a wide log-scale range,
    # LogLocator's default budget silently drops minor ticks entirely, and
    # thins major tick MARKS (not just labels) at every other decade.
    plt.gca().yaxis.set_minor_locator(
        ticker.LogLocator(base=10.0, subs=[2, 3, 4, 5, 6, 7, 8, 9], numticks=100))
    plt.gca().yaxis.set_minor_formatter(ticker.NullFormatter())
    plt.gca().yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=100))
    plt.gca().yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda val, pos: f"$10^{{{round(math.log10(val))}}}$"
                              if round(math.log10(val)) % 2 == 0 else ""))

    # Discrete per-y legend instead of a colorbar (see cross_section.py).
    y_handles = [Line2D([0], [0], color=colors[y], linestyle=y_linestyles[y], linewidth=2, label=f"$y={y:g}$")
                 for y in Y_TO_PLOT]
    plt.legend(handles=y_handles, loc="upper right", fontsize=15)

    plt.tight_layout()
    outname = f"../plots/cross_section_sum_{CHANNEL}_{NUCLEUS}.pdf"
    plt.savefig(outname, bbox_inches="tight")
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
