import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Reuses cross_section.py's data-loading (and picks up its "fancy" rcParams
# styling too, since that's set at import time). Plots the *total* D0 cross
# section (diffractive + exclusive added together) for a few rapidities.
from cross_section import load_results, NUCLEUS, CHANNEL, FRAG_TYPES

Y_TO_PLOT = [-1.0, 0.0, 1.0, 2.0, 3.0]


def main():
    frag = FRAG_TYPES[0]
    results_diff = load_results("diffractive", frag)
    results_excl = load_results("exclusive", frag)

    missing = [y for y in Y_TO_PLOT if y not in results_diff and y not in results_excl]
    if missing:
        sys.exit(f"No data found for rapidities {missing} -- run ../run_many_nucleus.sh first.")

    # Same ColorBrewer "Blues" scale as cross_section.py.
    blue_ramp = [
        "#deebf7", "#c6dbef", "#9ecae1", "#6baed6",
        "#4292c6", "#2171b5", "#08519c", "#08306b",
    ]
    colors = {}
    for i, y in enumerate(Y_TO_PLOT):
        step = round(i * (len(blue_ramp) - 1) / max(len(Y_TO_PLOT) - 1, 1))
        colors[y] = blue_ramp[step]

    plt.figure(figsize=(7.5, 6.5))
    for y in Y_TO_PLOT:
        diff_points = dict(results_diff.get(y, []))
        excl_points = dict(results_excl.get(y, []))
        pt_values = sorted(set(diff_points) | set(excl_points))
        totals = [diff_points.get(pt, 0.0) + excl_points.get(pt, 0.0) for pt in pt_values]
        plt.plot(pt_values, totals, color=colors[y], linestyle="-", linewidth=2)

    plt.yscale("log")
    plt.xlabel(r"$p_{D^0}$ [GeV]")
    plt.ylabel(r"$d\sigma/dy\,dp_{D^0}$ [mb/GeV]", labelpad=15)
    plt.title(f"$D^0$ photoproduction (diff. + excl.), {NUCLEUS}+{NUCLEUS} UPC ({CHANNEL})", pad=15)

    handles = [Line2D([0], [0], color=colors[y], linestyle="-", linewidth=2, label=f"y={y:g}")
               for y in Y_TO_PLOT]
    plt.legend(handles=handles, loc="upper right", fontsize=17)

    plt.tight_layout()
    outname = f"../plots/cross_section_sum_{CHANNEL}_{NUCLEUS}.pdf"
    plt.savefig(outname)
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
