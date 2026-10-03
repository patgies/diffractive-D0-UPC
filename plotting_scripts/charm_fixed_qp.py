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
    "axes.labelsize": 28,
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
    "xtick.color": "0.4",
    "ytick.color": "0.4",
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
    linestyle_cycle = ["-", "--", ":", "-."]
    linestyles = [linestyle_cycle[i % len(linestyle_cycle)] for i in range(len(y_values))]

    fig, ax = plt.subplots(figsize=(8, 6.5))

    data_by_y = {}
    for filename, y, linestyle in zip(files, y_values, linestyles):
        pt, excl, diff = read_data(filename)
        data_by_y[y] = (pt, excl, diff)
        ax.plot(pt, diff, color=DIFF_COLOR, linestyle=linestyle, lw=2.0)
        ax.plot(pt, excl, color=EXCL_COLOR, linestyle=linestyle, lw=2.0)

    ax.text(4.7, 4e-8, r"$\sim p_{\perp}^{-6}$", fontsize=22)
    ax.text(9.8, 5.5e-10, r"$\sim p_{\perp}^{-4}$", fontsize=22)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(bottom=1e-11)
    ax.set_xlim(0.2, 20)
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20))
    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20))
    ax.set_xlabel(r"$p_{\perp}$ [GeV]", labelpad=6)
    ax.set_ylabel(r"$d\sigma/dy\, d^2\mathbf{p}$ [mb/GeV$^2$]", labelpad=8)
    ax.set_title(r"$\gamma+A \rightarrow c+A$ at $z=0.5$", pad=15, fontsize=24)

    y_handles = [Line2D([0], [0], color='0.3', lw=3, linestyle=linestyles[i], label=f"$y={y:g}$")
                 for i, y in reversed(list(enumerate(y_values)))]
    style_handles = [Line2D([0], [0], color=EXCL_COLOR, lw=3, linestyle='-', label="Exclusive"),
                     Line2D([0], [0], color=DIFF_COLOR, lw=3, linestyle='-',
                            label=r"Diffractive$_{\mbox{\fontsize{13.5}{13.5}\selectfont TMD}}$")]
    y_legend = ax.legend(handles=y_handles, fontsize=22, loc='lower left',
                         bbox_to_anchor=(0.02, 0.0), frameon=False)
    ax.add_artist(y_legend)
    ax.legend(handles=style_handles, fontsize=22, loc='lower left',
              bbox_to_anchor=(0.33, 0.0), frameon=False)
    plt.tight_layout()
    outname = "../plots/charm_fixed_qp.pdf"
    plt.savefig(outname, dpi=150)
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
