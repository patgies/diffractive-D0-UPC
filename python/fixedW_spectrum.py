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

# Plots the fixed-W (no photon flux), q+=3p+ exclusive & diffractive spectra
# vs k_perp,c written by ../local_workflows/run_fixedW.sh (produced by the D0_fixedW
# program), one curve per rapidity.
#
# Normalization: diffractive already carries the full physical prefactor
# (see diffractiveCrossSection_fixedW in int_diffractive.cpp); exclusive
# still needs the prefactor.

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 20,
    "axes.titlesize": 20,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
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
    # differential, k_perp,c is fixed per point, not integrated); convert to
    # the radial dsigma/dK plotted here via the standard 2*pi*K Jacobian.
    jac = 2 * np.pi * K
    return K, jac * excl * PREFACTOR_EXCL, jac * diff


def main():
    pattern = "../files/D0_fixedW_y*.dat"
    files = sorted(glob.glob(pattern), key=read_y)
    if not files:
        sys.exit(f"No files found matching {pattern} -- run ../local_workflows/run_fixedW.sh first.")

    y_values = [read_y(f) for f in files]

    blue_ramp = ["#cde2fb", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
    colors = []
    for i in range(len(y_values)):
        step = round(i * (len(blue_ramp) - 1) / max(len(y_values) - 1, 1))
        colors.append(blue_ramp[step])

    fig, ax = plt.subplots(figsize=(8, 6.5))

    data_by_y = {}
    for filename, y, color in zip(files, y_values, colors):
        K, excl, diff = read_data(filename)
        data_by_y[y] = (K, excl, diff)
        ax.plot(K, diff, color=color, linestyle='-')
        ax.plot(K, excl, color=color, linestyle='--')

    # Power-law guide lines, anchored to the y=0 curves.
    K0, excl0, diff0 = data_by_y[0.0]
    K_anchor = K0[0]

    K_ref4 = np.array([K_anchor, 6.0])
    ax.plot(K_ref4, diff0[0] * (K_ref4/K_anchor)**(-4), color='black', linestyle=':', linewidth=1.5)
    ax.text(0.6, 4e-8, r"$\textrm{exclusive}\propto 1/p_{c\perp}^6$", fontsize=12,
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

    K_ref6 = np.array([K_anchor, 4.0])
    ax.plot(K_ref6, excl0[0] * (K_ref6/K_anchor)**(-6), color='black', linestyle=':', linewidth=1.5)
    ax.text(0.5, 3e-11, r"$\textrm{diffractive}\propto 1/p_{c\perp}^4$", fontsize=12,
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(bottom=1e-12)
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20))
    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    ax.set_xlabel(r"$p_{c\perp}$ [GeV]")
    ax.set_ylabel(r"$d\sigma/dy_c dp_{c\perp}$ [mb/GeV]", labelpad=15)
    ax.set_title(r"Charm photoproduction with proton target  (fixed $q^+=3p^+$)", pad=15)

    y_handles = [Line2D([0], [0], color=colors[i], linestyle='-', label=f"$y={y:g}$")
                 for i, y in reversed(list(enumerate(y_values)))]
    ax.legend(handles=y_handles, fontsize=13, loc='upper right',
              title="solid: diffractive\ndashed: exclusive", title_fontsize=12)
    plt.tight_layout()
    outname = "../plots/fixedW_spectrum.pdf"
    plt.savefig(outname, dpi=150)
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
