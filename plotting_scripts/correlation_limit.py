import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullFormatter

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 24,
    "axes.titlesize": 22,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
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

# Validity of the correlation limit: hard scale against the saturation scale, for a charm quark at fixed
# q+ = 2p+ (z = 1/2). Left, gluon DTMD: M_h = 2 m_T against the adjoint Q_s. Right, quark DTMD: p_perp
# (and m_c) against the fundamental Q_s. The bands span y = 0 to 2, with Q_s at the smallest x_po of each
# term (x < 1), where it is largest. Proton and Pb at b_d = 0.
# Inputs: output/saturation_scale/ (run_scripts/run_saturation_scale.sh).

QS_DIR = "../output/saturation_scale"
M_CHARM, SQRT_S = 1.5, 5360.0
XPO_MIN, XPO_FROZEN = 1e-6, 0.01   # as in def.hpp: Q_s is frozen above XPO_FROZEN


def qs_at(filename, x, column):
    """Q_s [GeV] at x from the saturation_scale output (column 1: fundamental, 2: adjoint)."""
    data = np.loadtxt(filename, comments='#')
    x = np.clip(x, XPO_MIN, XPO_FROZEN)
    return np.sqrt(np.exp(np.interp(np.log(x), np.log(data[:, 0]), np.log(data[:, column]))))


def xpo_min(p, y, M2):
    """Smallest x_po (x < 1) at q+ = 2p+: M^2/W^2 with W^2 = 2 m_T e^y sqrt(s)."""
    mt = np.hypot(p, M_CHARM)
    return M2 / (2 * mt * np.exp(y) * SQRT_S)


def main():
    TARGETS = [("proton", "qs_proton.dat", "#2166ac"), ("Pb", "qs_Pb_b0.dat", "#b2182b")]
    p = np.geomspace(0.2, 20, 200)
    mt = np.hypot(p, M_CHARM)

    fig, (left, right) = plt.subplots(1, 2, figsize=(13, 6), sharey=True)

    # Q_s bands (y = 0 to 2): gluon DTMD with M_qqbar^2 = 4 m_T^2, quark DTMD with M_qg^2 = 4 p^2 + 2 m^2
    for name, qs_file, color in TARGETS:
        qs_A = [qs_at(f"{QS_DIR}/{qs_file}", xpo_min(p, y, 4 * mt**2), 2) for y in (0.0, 2.0)]
        qs_F = [qs_at(f"{QS_DIR}/{qs_file}", xpo_min(p, y, 4 * p**2 + 2 * M_CHARM**2), 1) for y in (0.0, 2.0)]
        left.fill_between(p, qs_A[0], qs_A[1], color=color, alpha=0.35, lw=0)
        right.fill_between(p, qs_F[0], qs_F[1], color=color, alpha=0.35, lw=0)
        left.text(0.23, qs_A[1][0] * 1.08, rf"$Q_{{s,A}}$ {name}", fontsize=18, color=color)
        right.text(8, qs_F[1][-1] * 1.12, rf"$Q_s$ {name}", fontsize=18, color=color)

    left.plot(p, 2 * mt, color="k", lw=2.4)
    left.text(0.23, 2 * M_CHARM * 1.15, r"$M_h=2m_T$", fontsize=20)
    right.plot(p, p, color="k", lw=2.4)
    right.text(5, 6.6, r"$p_\perp$", fontsize=20)
    right.axhline(M_CHARM, color="0.3", lw=1.5, linestyle="--")
    right.text(0.23, M_CHARM * 1.08, r"$m_c$", fontsize=20, color="0.3")

    left.set_title("gluon DTMD", pad=12)
    right.set_title("quark DTMD", pad=12)
    for ax in (left, right):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(0.2, 20)
        ax.set_ylim(0.2, 50)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel(r"$p_\perp$ [GeV]", labelpad=6)
    left.set_ylabel("scale [GeV]", labelpad=6)

    plt.tight_layout()
    outname = "../plots/correlation_limit.pdf"
    plt.savefig(outname, dpi=150)
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
