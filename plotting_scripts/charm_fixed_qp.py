import glob
import math
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator

# Charm quark cross section at fixed q+ = 2p+ (z = 0.5) on a proton: no photon flux, no fragmentation.
# Reads the output of run_charm_fixed_qp.sh.
from D0 import alphas_run, alphae, mc, e_c, Nc, SIGMA0_MB   # also sets the plot style

plt.rcParams["axes.labelsize"] = 29

PREFACTOR_EXCL = alphae * Nc * e_c**2 * SIGMA0_MB / (2 * math.pi**2)


def prefactor_diff(pt):
    """Diffractive prefactor, with alpha_s at the scale m_T = sqrt(pt^2 + mc^2)."""
    alphas = alphas_run(np.sqrt(pt**2 + mc**2))
    return alphas * alphae * e_c**2 * (Nc**2 - 1) * SIGMA0_MB / (8 * math.pi**4)


def read_y(filename):
    with open(filename) as f:
        for line in f:
            if "y (fixed)" in line:
                return float(line.split("=")[-1])
    raise ValueError(f"no 'y (fixed)' header found in {filename}")


def read_data(filename):
    data = np.loadtxt(filename, comments='#')
    pt = data[:, 0]
    excl = data[:, 1]
    diff = data[:, 2]
    return pt, excl * PREFACTOR_EXCL, diff * prefactor_diff(pt)


def main():
    pattern = "../output/charm/charm_fixed_qp_y*.dat"
    files = sorted(glob.glob(pattern), key=read_y)
    if not files:
        sys.exit(f"No files found matching {pattern} -- run ../run_scripts/run_charm_fixed_qp.sh first.")

    Y_TO_PLOT = [0.0, 1.0, 2.0]   # rapidities shown (the files also have y = 0.5 and 1.5)
    files = [f for f in files if read_y(f) in Y_TO_PLOT]
    y_values = [read_y(f) for f in files]

    EXCL_COLOR, DIFF_COLOR = "#2166ac", "#b2182b"
    linestyle_cycle = ["-", "--", ":"]   # solid, dashed, dotted
    linestyles = [linestyle_cycle[i % len(linestyle_cycle)] for i in range(len(y_values))]

    fig, ax = plt.subplots(figsize=(8, 6.5))

    data_by_y = {}
    for filename, y, linestyle in zip(files, y_values, linestyles):
        pt, excl, diff = read_data(filename)
        data_by_y[y] = (pt, excl, diff)
        ax.plot(pt, diff, color=DIFF_COLOR, linestyle=linestyle, lw=2.0)
        ax.plot(pt, excl, color=EXCL_COLOR, linestyle=linestyle, lw=2.0)

    # guide lines with their labels: p^-8 under the exclusive curves, p^-4 above the diffractive ones
    # (between the two sets of curves there is no room for a label)
    x8 = np.geomspace(4.2, 7.8, 40)
    ax.plot(x8, 7e-9 * (x8 / 4) ** -8, color="0.5", lw=1.8, zorder=1)
    x4 = np.geomspace(4, 9.5, 20)
    ax.plot(x4, 1.3e-6 * (x4 / 5) ** -4, color="0.5", lw=1.8, zorder=1)
    ax.text(2.95, 3.5e-10, r"$\sim p_{\perp}^{-8}$", fontsize=23)
    ax.text(5.2, 1.6e-6, r"$\sim p_{\perp}^{-4}$", fontsize=23)

    # log-log, so that the power laws are straight lines
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e-11, 2e-4)
    ax.set_xlim(0.2, 10)
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20))
    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20))
    ax.set_xlabel(r"$p_{\perp}$ [GeV]", labelpad=6)
    ax.set_ylabel(r"$d\sigma/dy\, d^2\mathbf{p}$ [mb/GeV$^2$]", labelpad=4, fontsize=27)
    ax.set_title(r"$\gamma+p \rightarrow c+p$ at $z=0.5$", pad=15, fontsize=28)

    y_handles = [Line2D([0], [0], color='black', lw=3, linestyle=linestyles[i], label=f"$y={y:g}$")
                 for i, y in reversed(list(enumerate(y_values)))]
    style_handles = [Line2D([0], [0], color=EXCL_COLOR, lw=3, linestyle='-', label="Exclusive"),
                     Line2D([0], [0], color=DIFF_COLOR, lw=3, linestyle='-',
                            label=r"Diffractive$_{\mbox{\fontsize{14}{14}\selectfont TMD}}$")]
    y_legend = ax.legend(handles=y_handles, fontsize=21, loc='lower left',
                         bbox_to_anchor=(0.02, 0.0), frameon=False)
    ax.add_artist(y_legend)
    ax.legend(handles=style_handles, fontsize=23, loc='lower left', bbox_to_anchor=(0.02, 0.27), frameon=False)
    plt.tight_layout()
    outname = "../plots/charm_fixed_qp.pdf"
    plt.savefig(outname, dpi=150)
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
