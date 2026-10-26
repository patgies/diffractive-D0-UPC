import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator

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
    "xtick.color": "gray",
    "ytick.color": "gray",
    "xtick.labelcolor": "black",
    "ytick.labelcolor": "black",
    "xtick.major.pad": 8,
    "ytick.major.pad": 4,
})
plt.rcParams["text.latex.preamble"] = r"\usepackage{amssymb}"   # \mathbb{P}

# Saturation scale Q_s = sqrt(2)/r_s of the proton and of Pb vs x_po, for the fundamental dipole N
# (dashed) and the adjoint dipole N_A = 2N - N^2 of the diffractive part (solid), to compare with
# the p_perp of the predictions (the diffractive TMD needs p_perp >> Q_s). Pb: centre of the nucleus
# (b_d = 0). Above x_po = 0.01 the dipole is frozen.
# The vertical lines are the smallest x_po of the Pb+Pb diffractive part at y = 0 and y = 2 (x < 1 at
# z_h = 1, the smallest pD0 and q+ = qpmax); in the gray region x_po is not reached at |y| <= 2.
# Inputs: output/saturation_scale/, made by run_scripts/run_saturation_scale.sh.

DATA_DIR = "../output/saturation_scale"
M_CHARM, SQRT_S, QP_MAX, PD0_MIN = 1.5, 5360.0, 800.0, 0.2   # as in setup.cpp and config.sh


def xpo_min(y):
    """Smallest x_po with x = M_qq^2/(x_po W^2) < 1: z_h = 1, pD0 = PD0_MIN and q+ = QP_MAX."""
    mt = np.hypot(PD0_MIN, M_CHARM)
    pp = mt / np.sqrt(2) * np.exp(y)
    return mt**2 * QP_MAX / (pp * (QP_MAX - pp) * np.sqrt(2) * SQRT_S)


def read_qs(filename):
    """x, Q_s fundamental, Q_s adjoint [GeV]."""
    data = np.loadtxt(filename, comments='#')
    return data[:, 0], np.sqrt(data[:, 1]), np.sqrt(data[:, 2])


def main():
    P_COLOR, PB_COLOR = "#2166ac", "#b2182b"

    fig, ax = plt.subplots(figsize=(8, 6.5))

    x_lo = 1e-5
    for target, color in [("proton", P_COLOR), ("Pb_b0", PB_COLOR)]:
        x, qs_f, qs_a = read_qs(f"{DATA_DIR}/qs_{target}.dat")
        ax.plot(x, qs_a, color=color, lw=2.0)
        ax.plot(x, qs_f, color=color, lw=2.0, linestyle="--")

    ax.axvspan(x_lo, xpo_min(2.0), color="0.85", lw=0, zorder=0)
    for y in [0.0, 2.0]:
        ax.axvline(xpo_min(y), color="0.4", lw=1.5, linestyle=":", zorder=1)
        ax.text(xpo_min(y) * 1.12, 0.15, f"$y={y:g}$", fontsize=20, color="0.3")

    ax.set_xscale('log')
    ax.set_xlim(x_lo, 1e-2)
    ax.set_ylim(0, 4.5)
    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20))
    ax.set_xlabel(r"$x_{\mathbb{P}}$", labelpad=6)
    ax.set_ylabel(r"$Q_s$ [GeV]", labelpad=6)
    ax.set_title(r"$N(r_s) = 1 - e^{-1/2}$, $\; Q_s^2 = 2/r_s^2$", pad=15, fontsize=24)

    target_handles = [Line2D([0], [0], color=PB_COLOR, lw=3, label=r"Pb, $b_d = 0$"),
                      Line2D([0], [0], color=P_COLOR, lw=3, label="proton")]
    style_handles = [Line2D([0], [0], color='black', lw=3, linestyle='-', label=r"adjoint $N_A$"),
                     Line2D([0], [0], color='black', lw=3, linestyle='--', label=r"fundamental $N$")]
    target_legend = ax.legend(handles=target_handles, fontsize=21, loc='upper right', frameon=False)
    ax.add_artist(target_legend)
    ax.legend(handles=style_handles, fontsize=21, loc='upper right', bbox_to_anchor=(1.0, 0.76), frameon=False)
    plt.tight_layout()
    outname = "../plots/saturation_scale.pdf"
    plt.savefig(outname, dpi=150)
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
