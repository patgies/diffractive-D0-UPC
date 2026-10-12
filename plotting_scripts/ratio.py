import glob
import math
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator, RegularGridInterpolator

from D0 import (
    read_rapidity, read_data_file, group_by_pt, integrate_over_b_d,
    load_results, _summed_results, alphae, e_c, Nc, GEVSQR_TO_MB, NUCLEUS, CHANNEL, CENTRAL_DIR,
)

# Contour plot of (diffractive+exclusive)/inclusive D0 in the (y, pT) plane, for the fragmentation
# function FRAG_TYPE (default BCFY). Inclusive files: ../input/inclusive/, from inclusive-D0-UPC.

INCLUSIVE_FACTOR_A = alphae * e_c**2 * Nc / (2 * math.pi) ** 4   # prefactor of inclusive-D0-UPC

PT_MAX_PLOT = 8.0


def load_inclusive_results(frag="HymnD"):
    """Same file format as load_results, but with inclusive-D0-UPC's own normalization."""
    pattern = f"../input/inclusive/D0_incl_{frag}_{CHANNEL}_{NUCLEUS}_y*.dat"
    results = {}
    for filename in sorted(glob.glob(pattern)):
        y = read_rapidity(filename)
        b_d_list, pt_list, dsigma_list = read_data_file(filename)
        pt_groups = group_by_pt(b_d_list, pt_list, dsigma_list)

        points = []
        for pt, pairs in sorted(pt_groups.items()):
            b_d_integral = integrate_over_b_d(pairs)
            cross_section = (2*math.pi) * b_d_integral * INCLUSIVE_FACTOR_A \
                            * (2*math.pi) * pt * GEVSQR_TO_MB
            points.append((pt, cross_section))
        if points:
            results[y] = points
    return results


def compute_ratio_grid(process_results, incl_results, y_values, pt_grid):
    """Returns a (len(pt_grid), len(y_values)) array of process/inclusive, with the
    inclusive spectrum interpolated to pt_grid."""
    ratio_grid = np.full((len(pt_grid), len(y_values)), np.nan)
    for iy, y in enumerate(y_values):
        proc_points = dict(process_results[y])
        incl_points = dict(incl_results[y])
        incl_pt = np.array(sorted(incl_points))
        # smooth (cubic) interpolation of ln(dsigma/d^2p): without the factor pt the spectrum has no
        # maximum, and a linear interpolation gives a sawtooth in the ratio (the two pt grids differ)
        incl_spline = PchipInterpolator(incl_pt, np.log(np.array([incl_points[pt] for pt in incl_pt]) / incl_pt))

        for ipt, pt in enumerate(pt_grid):
            if pt not in proc_points:
                continue
            incl_xs = pt * math.exp(float(incl_spline(pt)))
            ratio_grid[ipt, iy] = proc_points[pt] / incl_xs

    return ratio_grid


def label_positions(y_values, pt_grid, ratio_grid, levels):
    """One label position per contour level, where the level crosses the diagonal from the
    upper left to the lower right of the panel (so the labels do not overlap)."""
    interpolate = RegularGridInterpolator((np.array(pt_grid), np.array(y_values)), ratio_grid)
    t = np.linspace(0.0, 1.0, 400)
    y_line = (y_values[0] + 0.25) + t * ((y_values[-1] - 0.9) - (y_values[0] + 0.25))
    pt_line = (pt_grid[-1] - 0.4) + t * ((pt_grid[0] + 0.35) - (pt_grid[-1] - 0.4))
    values = interpolate(np.column_stack([pt_line, y_line]))
    positions = []
    for level in levels:
        crossings = np.where(np.diff(np.sign(values - level)) != 0)[0]
        if len(crossings):
            positions.append((y_line[crossings[0]], pt_line[crossings[0]]))
    return positions


def plot_ratio_panel(ax, y_values, pt_grid, ratio_grid, label):
    Y, PT = np.meshgrid(y_values, pt_grid)

    contourf = ax.contourf(Y, PT, ratio_grid, levels=20, cmap="Blues")
    line_levels = [0.04, 0.05, 0.06, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30, 0.35]
    contour_lines = ax.contour(Y, PT, ratio_grid, levels=line_levels, colors="navy",
                                linestyles="-", linewidths=0.6, alpha=0.7)
    # labels in black, bold and horizontal
    contour_labels = ax.clabel(contour_lines, inline=True, fontsize=14, colors="black",
                               fmt=lambda v: r"\textbf{%.2f}" % v,
                               manual=label_positions(y_values, pt_grid, ratio_grid, line_levels))
    for contour_label in contour_labels:
        contour_label.set_rotation(0)
        contour_label.set_alpha(1.0)
        contour_label.set_color("black")

    colorbar = plt.colorbar(contourf, ax=ax)
    ticks = colorbar.get_ticks()
    colorbar.set_ticks(ticks)
    colorbar.ax.set_yticklabels(["0" if abs(tick) < 1e-9 else f"{tick:.2f}" for tick in ticks])   # 0, not 0.00
    colorbar.ax.tick_params(which="both", color="black", labelsize=25)
    ax.tick_params(which="both", color="black")

    ax.set_xlabel(r"$y$", fontsize=28)
    ax.set_ylabel(r"$p_{D\perp}$ [GeV]", fontsize=28, labelpad=10)
    ax.set_title(label, pad=12, fontsize=24, linespacing=1.5)
    ax.tick_params(axis="both", labelsize=25)
    ax.tick_params(axis="y", pad=9)
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
    frag = os.environ.get("FRAG_TYPE", "BCFY")
    process = "sum"

    proc_results = load_process_results(process, frag)
    incl_results = load_inclusive_results(frag)

    y_values = sorted(set(proc_results) & set(incl_results))
    if not y_values:
        sys.exit(f"No overlapping rapidities across {process}/inclusive data for {frag} -- "
                  f"check {CENTRAL_DIR} and ../input/inclusive/.")
    y_range = (y_values[0], y_values[-1], len(y_values))

    pt_grid = sorted({pt for pt, _ in proc_results[y_values[0]] if pt <= PT_MAX_PLOT})
    ratio = compute_ratio_grid(proc_results, incl_results, y_values, pt_grid)

    fig, ax = plt.subplots(figsize=(8, 6.5))
    # "Diffractive" = exclusive + diffractive TMD, as in the other figures
    title = r"$\sigma^{\rm diffractive}/\sigma^{\rm inclusive}$"
    plot_ratio_panel(ax, y_values, pt_grid, ratio, title)

    plt.tight_layout()
    outname = f"../plots/ratio_{CHANNEL}_{NUCLEUS}{NUCLEUS}_{frag}.pdf"
    plt.savefig(outname)
    print(f"Saved: {outname}")
    print(f"y range used: {y_range[0]:g} to {y_range[1]:g} ({y_range[2]} points)")


if __name__ == "__main__":
    main()
