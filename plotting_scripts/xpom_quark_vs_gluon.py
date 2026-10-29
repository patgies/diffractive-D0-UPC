import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Gluon vs quark DTMD x_po distribution at fixed q+ = 2p+ (z = 0.5), p = 8 GeV, y = 0: no photon flux,
# no fragmentation, proton target as in charm_fixed_qp_xpom.py. Data from
#   DTMD=gluon|quark ./build/bin/charm_fixed_qp_xpom 0 8 30 1e5
from charm_fixed_qp import prefactor_diff   # also sets the plot style

plt.rcParams["text.latex.preamble"] = r"\usepackage{amssymb}"   # \mathbb{P}
P = 8.0


def load(dtmd):
    """Returns (x_po, dsigma/dy d2p dln x_po) above the smallest x_po, and that x_po."""
    filename = f"../output/charm/xpom/charm_fixed_qp_xpom_{dtmd}_y00_p80.dat"
    with open(filename) as f:
        xpo_min = next(float(line.split("=")[1].split()[0]) for line in f if line.startswith("# x_po (excl)"))
    data = np.loadtxt(filename, comments="#")
    keep = data[:, 1] > 0
    return data[keep, 0], data[keep, 1] * prefactor_diff(P), xpo_min


def main():
    xg, g, xpo_min = load("gluon")
    xq, q, _ = load("quark")

    plt.figure(figsize=(8, 6.5))
    plt.plot(xg, g, color="#b2182b", lw=2.5, label="gluon DTMD")
    plt.plot(xq, q, color="#4dac26", lw=2.5, label="quark DTMD")
    plt.axvline(xpo_min, color="0.5", ls=":", lw=2, label=r"$x_{\mathbb{P}}^{\min}=M^2/W^2$")

    plt.xscale("log")
    plt.yscale("log")
    plt.ylim(3e-12, 1e-7)
    plt.xlabel(r"$x_{\mathbb{P}}$", labelpad=10)
    plt.ylabel(r"$d\sigma/dy\,d^2\mathbf{p}\,d\ln(x_{\mathbb{P}})$ [mb/GeV$^2$]", labelpad=15, fontsize=27)
    plt.title(rf"$\gamma+p$, $p_\perp={P:g}$ GeV, $y=0$, $z=0.5$", pad=15, fontsize=28)
    plt.legend(loc="center right", fontsize=23, frameon=False)

    plt.tight_layout()
    outname = "../plots/xpom_quark_vs_gluon.pdf"
    plt.savefig(outname, bbox_inches="tight")
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
