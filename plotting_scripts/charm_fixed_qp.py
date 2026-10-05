import glob
import math
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))   # alphas_running.py
from alphas_running import alphas_run


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 29,
    "axes.titlesize": 20,
    "xtick.labelsize": 25,
    "ytick.labelsize": 25,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 10,
    "ytick.major.size": 10,
    "xtick.minor.size": 5,
    "ytick.minor.size": 5,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.color": "gray",
    "ytick.color": "gray",
    "xtick.labelcolor": "black",
    "ytick.labelcolor": "black",
    "xtick.major.pad": 8,
    "ytick.major.pad": 4,
})

alphae = 1/137
mc     = 1.5
e_c    = 2/3
Nc     = 3
sigma0 = 16.36

PREFACTOR_EXCL = alphae * Nc * e_c**2 * sigma0 / (2 * math.pi**2)


def prefactor_diff(pt):
    """Diffractive prefactor, with alpha_s at the scale m_T = sqrt(pt^2 + mc^2)."""
    alphas = alphas_run(np.sqrt(pt**2 + mc**2))
    return alphas * alphae * e_c**2 * (Nc**2 - 1) * sigma0 / (8 * math.pi**4)


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

    y_values = [read_y(f) for f in files]

  
    EXCL_COLOR, DIFF_COLOR = "#2166ac", "#b2182b"
    linestyle_cycle = ["-", (0, (8, 2)), "--", "-.", (0, (5, 1.5, 1, 1.5, 1, 1.5)), ":"]   # as in RpA_fixed_qp.py
    linestyles = [linestyle_cycle[i % len(linestyle_cycle)] for i in range(len(y_values))]

    fig, ax = plt.subplots(figsize=(8, 6.5))

    data_by_y = {}
    for filename, y, linestyle in zip(files, y_values, linestyles):
        pt, excl, diff = read_data(filename)
        data_by_y[y] = (pt, excl, diff)
        ax.plot(pt, diff, color=DIFF_COLOR, linestyle=linestyle, lw=2.0)
        ax.plot(pt, excl, color=EXCL_COLOR, linestyle=linestyle, lw=2.0)

    x6 = np.geomspace(4, 8.5, 40)
    ax.plot(x6, 6e-8 * (x6 / 4.6) ** -6, color="0.5", lw=1.8, zorder=1)
    x4 = np.geomspace(8.5, 20, 20)
    ax.plot(x4, 4e-10 * (x4 / 10) ** -4, color="0.5", lw=1.8, zorder=1)
    ax.text(5.1, 4.5e-8, r"$\sim p_{\perp}^{-6}$", fontsize=23)
    ax.text(10.3, 5.5e-10, r"$\sim p_{\perp}^{-4}$", fontsize=23)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e-11, 1e-4)
    ax.set_xlim(0.9, 20)
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20))
    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20))
    ax.set_xlabel(r"$p_{\perp}$ [GeV]", labelpad=6)
    ax.set_ylabel(r"$d\sigma/dy\, d^2\mathbf{p}$ [mb/GeV$^2$]", labelpad=8)
    ax.set_title(r"$\gamma+A \rightarrow c+A$ at $z=0.5$", pad=15, fontsize=25)

    y_handles = [Line2D([0], [0], color='black', lw=3, linestyle=linestyles[i], label=f"$y={y:g}$")
                 for i, y in reversed(list(enumerate(y_values)))]
    style_handles = [Line2D([0], [0], color=EXCL_COLOR, lw=3, linestyle='-', label="Exclusive"),
                     Line2D([0], [0], color=DIFF_COLOR, lw=3, linestyle='-',
                            label=r"Diffractive$_{\mbox{\fontsize{14}{14}\selectfont TMD}}$")]
    y_legend = ax.legend(handles=y_handles, fontsize=23, loc='lower left',
                         bbox_to_anchor=(0.02, 0.0), frameon=False)
    ax.add_artist(y_legend)
    ax.legend(handles=style_handles, fontsize=23, loc='upper right', bbox_to_anchor=(1.0, 0.96), frameon=False)
    plt.tight_layout()
    outname = "../plots/charm_fixed_qp.pdf"
    plt.savefig(outname, dpi=150)
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
