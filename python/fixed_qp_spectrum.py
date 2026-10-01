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

# Plots the fixed-q+ (no photon flux), q+=2p+ exclusive & diffractive spectra
# vs k_perp,c written by ../local_workflows/run_fixed_qp.sh (produced by the D0_fixed_qp
# program), one curve per rapidity.
#
# Normalization: diffractive already carries the full physical prefactor
# (see diffractiveCrossSection_fixed_qp in int_diffractive.cpp); exclusive
# still needs the prefactor.

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 22,
    "axes.titlesize": 20,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 8,
    "ytick.major.size": 8,
    "xtick.minor.size": 4,
    "ytick.minor.size": 4,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    # same tick styling as flux_comparison.py: gray ticks, black numbers
    "xtick.color": "0.4",
    "ytick.color": "0.4",
    "xtick.labelcolor": "black",
    "ytick.labelcolor": "black",
    "xtick.major.pad": 8,
    "ytick.major.pad": 4,
})

alphae = 1/137
mc     = 1.5          # charm mass in GeV
e_c    = 2/3
Nc     = 3
sigma0 = 16.36   # dipole normalization, same convention as python/*.py

# Missing prefactor for the exclusive column (diffractive already has its
# full prefactor baked in on the C++ side, see int_diffractive.cpp).
PREFACTOR_EXCL = alphae * Nc * e_c**2 * sigma0 / (2 * math.pi**2)


def read_y(filename):
    with open(filename) as f:
        for line in f:
            if "y (fixed)" in line:
                return float(line.split("=")[-1])
    raise ValueError(f"no 'y (fixed)' header found in {filename}")


def read_data(filename):
    data = np.loadtxt(filename, comments='#')
    K = data[:, 0]
    excl = data[:, 1]
    diff = data[:, 2]
    # Both columns come out as dsigma/d^2K (2D transverse-momentum
    # differential, k_perp,c is fixed per point, not integrated), which is
    # what is plotted -- the large-K power laws 1/K^6 (exclusive) and 1/K^4
    # (diffractive) refer to this distribution.
    return K, excl * PREFACTOR_EXCL, diff


def main():
    pattern = "../files/D0_fixed_qp_y*.dat"
    files = sorted(glob.glob(pattern), key=read_y)
    if not files:
        sys.exit(f"No files found matching {pattern} -- run ../local_workflows/run_fixed_qp.sh first.")

    y_values = [read_y(f) for f in files]

    # Same red / green / blue as the channels in flux_comparison.py.
    palette = ["#d62728", "#2563b8", "#1a7f4b"]
    colors = [palette[i % len(palette)] for i in range(len(y_values))]

    fig, ax = plt.subplots(figsize=(8, 6.5))

    data_by_y = {}
    for filename, y, color in zip(files, y_values, colors):
        K, excl, diff = read_data(filename)
        data_by_y[y] = (K, excl, diff)
        ax.plot(K, diff, color=color, linestyle='--')
        ax.plot(K, excl, color=color, linestyle='-')

    # Power-law behaviour of the tails, labelled next to them.
    ax.text(4.7, 4e-8, r"$\sim p_{\perp}^{-6}$", fontsize=20)
    ax.text(11.2, 1.8e-10, r"$\sim p_{\perp}^{-4}$", fontsize=20)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(bottom=1e-11)
    ax.set_xlim(0.2, 20)   # the p_perp grid of run_fixed_qp.sh
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20))
    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20))
    ax.set_xlabel(r"$p_{\perp}$ [GeV]", labelpad=6)
    ax.set_ylabel(r"$d\sigma/dy\, d^2\mathbf{p}$ [mb/GeV$^2$]", labelpad=8)
    ax.set_title(r"Charm photoproduction at fixed $z=p^+/q^+=1/2$", pad=15, fontsize=21)

    # Two frameless legends, as in flux_comparison.py: color = rapidity,
    # line style (black) = exclusive / diffractive.
    y_handles = [Line2D([0], [0], color=colors[i], lw=2.5, linestyle='-', label=f"$y={y:g}$")
                 for i, y in reversed(list(enumerate(y_values)))]
    style_handles = [Line2D([0], [0], color='k', lw=2, linestyle='-', label="Exclusive"),
                     Line2D([0], [0], color='k', lw=2, linestyle='--', label="Diffractive")]
    y_legend = ax.legend(handles=y_handles, fontsize=20, loc='lower left',
                         bbox_to_anchor=(0.02, 0.0), frameon=False)
    ax.add_artist(y_legend)
    ax.legend(handles=style_handles, fontsize=20, loc='lower left',
              bbox_to_anchor=(0.3, 0.0), frameon=False)
    plt.tight_layout()
    outname = "../plots/fixed_qp_spectrum.pdf"
    plt.savefig(outname, dpi=150)
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
