import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullFormatter

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 28,
    "axes.titlesize": 26,
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

# Correlation limit check: hard scale vs Q_s for charm at q+ = 2p+. Left: gluon DTMD, 2 m_T vs adjoint Q_s.
# Right: quark DTMD, p vs fundamental Q_s. Bands: y = 0 to 3. Proton and Pb at b_d = 0.
# Input: output/saturation_scale/ (run_saturation_scale.sh).

QS_DIR = "../output/saturation_scale"
M_CHARM, SQRT_S = 1.5, 5360.0
XPO_MIN, XPO_FROZEN = 1e-6, 0.01   # Q_s frozen above XPO_FROZEN (def.hpp)


def qs_at(filename, x, column):
    """Q_s [GeV] at x (column 1: fundamental, 2: adjoint)."""
    data = np.loadtxt(filename, comments='#')
    x = np.clip(x, XPO_MIN, XPO_FROZEN)
    return np.sqrt(np.exp(np.interp(np.log(x), np.log(data[:, 0]), np.log(data[:, column]))))


def xpo_min(p, y, M2):
    """Smallest x_po at q+ = 2p+."""
    mt = np.hypot(p, M_CHARM)
    return M2 / (2 * mt * np.exp(y) * SQRT_S)


def main():
    TARGETS = [("proton", "qs_proton.dat", "#2166ac"), ("Pb", "qs_Pb_b0.dat", "#b2182b")]
    p = np.geomspace(0.2, 20, 200)
    mt = np.hypot(p, M_CHARM)

    fig, (left, right) = plt.subplots(1, 2, figsize=(13, 6), sharey=True)

    left.plot(p, 2 * mt, color="k", lw=2.4, zorder=3, label=r"$2\sqrt{p_\perp^2+m_c^2}$")
    right.plot(p, p, color="k", lw=2.4, zorder=3, label=r"$p_\perp$")
    right.axhline(M_CHARM, color="0.3", lw=1.5, linestyle="--", zorder=3, label=r"$m_c$")

    # Q_s bands: gluon with M^2 = 4 m_T^2, quark with M^2 = 4 p^2 + 2 m^2
    for name, qs_file, color in reversed(TARGETS):
        qs_A = [qs_at(f"{QS_DIR}/{qs_file}", xpo_min(p, y, 4 * mt**2), 2) for y in (0.0, 3.0)]
        qs_F = [qs_at(f"{QS_DIR}/{qs_file}", xpo_min(p, y, 4 * p**2 + 2 * M_CHARM**2), 1) for y in (0.0, 3.0)]
        left.fill_between(p, qs_A[0], qs_A[1], color=color, alpha=0.45, lw=0, label=rf"$Q_{{s,A}}$ {name}")
        right.fill_between(p, qs_F[0], qs_F[1], color=color, alpha=0.45, lw=0, label=rf"$Q_s$ {name}")

    for ax, title in ((left, "Gluon DTMD"), (right, "Quark DTMD")):
        ax.set_title(title, pad=15)
        ax.legend(loc="upper left", fontsize=22, frameon=False)
    for ax in (left, right):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(0.2, 20)
        ax.set_ylim(0.2, 50)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel(r"$p_\perp$ [GeV]", labelpad=6)
    left.set_ylabel("scale [GeV]", labelpad=14)

    plt.tight_layout()
    outname = "../plots/correlation_limit.pdf"
    plt.savefig(outname, dpi=150)
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
