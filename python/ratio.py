import glob
import math
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from D0 import (
    read_rapidity, read_data_file, group_by_pt, integrate_over_b,
    load_results, _summed_results, GEVSQR_TO_MB, NUCLEUS, CHANNEL, CENTRAL_DIR,
)

# Contour plot of (diffractive+exclusive)/inclusive D0 in the (y, pT) plane.
# Inclusive files: ../inputs/inclusive/.

alpha_em = 1/137
e_charm_squared = 4/9
Nc = 3
INCLUSIVE_FACTOR_A = alpha_em * e_charm_squared * Nc / (2 * math.pi) ** 4

PT_MAX_PLOT = 2.0


def load_inclusive_results(frag="HymnD"):
    """Same file format as load_results, but with inclusive-D0-UPC's own normalization."""
    pattern = f"../inputs/inclusive/D0_incl_{frag}_{CHANNEL}_{NUCLEUS}_y*.dat"
    results = {}
    for filename in sorted(glob.glob(pattern)):
        y = read_rapidity(filename)
        b_list, pt_list, dsigma_list = read_data_file(filename)
        pt_groups = group_by_pt(b_list, pt_list, dsigma_list)

        points = []
        for pt, pairs in sorted(pt_groups.items()):
            b_integral = integrate_over_b(pairs)
            cross_section = (2*math.pi) * b_integral * INCLUSIVE_FACTOR_A \
                            * (2*math.pi) * pt * GEVSQR_TO_MB
            points.append((pt, cross_section))
        if points:
            results[y] = points
    return results


def compute_ratio_grid(process_results, incl_results, y_values, pt_grid):
    """Returns a (len(pt_grid), len(y_values)) array of process/inclusive, with the
    inclusive spectrum interpolated (in log-log) to pt_grid."""
    ratio_grid = np.full((len(pt_grid), len(y_values)), np.nan)
    for iy, y in enumerate(y_values):
        proc_points = dict(process_results[y])
        incl_points = dict(incl_results[y])
        incl_pt = np.array(sorted(incl_points))
        incl_log_xs = np.log(np.array([incl_points[pt] for pt in incl_pt]))

        for ipt, pt in enumerate(pt_grid):
            if pt not in proc_points:
                continue
            log_incl_xs = np.interp(pt, incl_pt, incl_log_xs)
            incl_xs = math.exp(log_incl_xs)
            ratio_grid[ipt, iy] = proc_points[pt] / incl_xs

    return ratio_grid


def plot_ratio_panel(ax, y_values, pt_grid, ratio_grid, label):
    Y, PT = np.meshgrid(y_values, pt_grid)

    contourf = ax.contourf(Y, PT, ratio_grid, levels=20, cmap="Blues")
    contour_lines = ax.contour(Y, PT, ratio_grid, levels=10, colors="navy",
                                linestyles=":", linewidths=0.8)
    ax.clabel(contour_lines, inline=True, fontsize=13, fmt="%.2f")

    plt.colorbar(contourf, ax=ax)

    ax.set_xlabel(r"$y$", fontsize=24)
    ax.set_ylabel(r"$p_{D^0\perp}$ [GeV]", fontsize=24, labelpad=15)
    ax.set_title(label, pad=12, fontsize=24)
    ax.tick_params(axis="both", labelsize=20)
    plt.setp(ax.get_xticklabels(), fontweight="normal")
    plt.setp(ax.get_yticklabels(), fontweight="normal")
    ax.set_xlim(-3, 3)
    ax.set_xticks([-3, -2, -1, 0, 1, 2, 3])


def load_process_results(process, frag):
    """process is "diffractive"/"exclusive" (-> load_results) or "sum"
    (-> diffractive+exclusive per (y,pT), via _summed_results)."""
    if process == "sum":
        return _summed_results(CENTRAL_DIR, frag)
    return load_results(process, frag)


def main():
    frag = "BCFY"
    process = "sum"

    proc_results = load_process_results(process, frag)
    incl_results = load_inclusive_results(frag)

    y_values = sorted(set(proc_results) & set(incl_results))
    if not y_values:
        sys.exit(f"No overlapping rapidities across {process}/inclusive data for {frag} -- "
                  f"check {CENTRAL_DIR} and ../inputs/inclusive/.")
    y_range = (y_values[0], y_values[-1], len(y_values))

    pt_grid = sorted({pt for pt, _ in proc_results[y_values[0]] if pt <= PT_MAX_PLOT})
    ratio = compute_ratio_grid(proc_results, incl_results, y_values, pt_grid)

    fig, ax = plt.subplots(figsize=(8, 6.5))
    plot_ratio_panel(ax, y_values, pt_grid, ratio, f"(diffractive + exclusive) / inclusive ({frag})")

    plt.tight_layout()
    outname = f"../plots/ratio_{CHANNEL}_{NUCLEUS}{NUCLEUS}.pdf"
    plt.savefig(outname)
    print(f"Saved: {outname}")
    print(f"y range used: {y_range[0]:g} to {y_range[1]:g} ({y_range[2]} points)")


if __name__ == "__main__":
    main()
