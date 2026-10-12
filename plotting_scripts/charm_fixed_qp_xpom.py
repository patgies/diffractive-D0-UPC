import glob
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Diffractive charm x_po distribution at fixed q+ = 2p+ (z = 0.5): no photon flux, no fragmentation,
# proton target as in charm_fixed_qp.py. Reads the output of run_charm_fixed_qp_xpom.sh.
from charm_fixed_qp import prefactor_diff

plt.rcParams["text.latex.preamble"] = r"\usepackage{amssymb}"   # \mathbb{P}


def read_header(filename):
    """Returns (y, p) from the '# y (fixed)' and '# p (fixed)' lines."""
    values = {}
    with open(filename) as f:
        for line in f:
            for key in ("y", "p"):
                if line.startswith(f"# {key} (fixed)"):
                    values[key] = float(line.split("=")[-1])
    return values["y"], values["p"]


def main():
    pattern = "../output/charm/xpom/charm_fixed_qp_xpom_y*_p*.dat"
    files = sorted(glob.glob(pattern), key=lambda f: read_header(f)[0])
    if not files:
        sys.exit(f"No files found matching {pattern} -- run ../run_scripts/run_charm_fixed_qp_xpom.sh first.")

    Y_TO_PLOT = [0.0, 0.5, 1.0, 1.5, 2.0]
    files = [f for f in files if read_header(f)[0] in Y_TO_PLOT]

    y_colors = {0.0: "#2166ac", 0.5: "#67a9cf", 1.0: "#ef8a62", 1.5: "#d6604d", 2.0: "#b2182b"}
    linestyle_cycle = ["-", (0, (8, 2)), "--", "-.", ":"]

    plt.figure(figsize=(8, 7))
    handles = []
    p_values = set()
    for i, filename in enumerate(files):
        y, p = read_header(filename)
        p_values.add(p)
        data = np.loadtxt(filename, comments="#")
        xpo, dsigma_dlnxpo = data[:, 0], data[:, 1] * prefactor_diff(p)
        keep = dsigma_dlnxpo > 0   # zero below the exclusive x_po
        color, linestyle = y_colors.get(y, "0.3"), linestyle_cycle[i % len(linestyle_cycle)]
        plt.plot(np.log(xpo[keep]), dsigma_dlnxpo[keep], color=color, linestyle=linestyle, lw=2.0)
        handles.append(Line2D([0], [0], color=color, linestyle=linestyle, lw=3, label=f"$y={y:g}$"))

    plt.yscale("log")
    plt.ylim(3e-9, 1e-4)
    plt.xlim(-9, np.log(0.1))
    plt.xlabel(r"$\ln(x_{\mathbb{P}})$", labelpad=14)
    plt.ylabel(r"$d\sigma/dy\,d^2\mathbf{p}\,d\ln(x_{\mathbb{P}})$ [mb/GeV$^2$]", labelpad=15, fontsize=27)

    ax = plt.gca()
    ax.tick_params(axis='x', which='both', top=False)
    ax.tick_params(axis='y', which='both', right=True)
    secax = ax.secondary_xaxis('top', functions=(np.exp, np.log))
    secax.set_xscale('log')
    secax.set_xlabel(r"$x_{\mathbb{P}}$", labelpad=12)

    p_label = ", ".join(f"{p:g}" for p in sorted(p_values))
    ax.text(0.08, 0.91, r"Diffractive$_{\mbox{\fontsize{14}{14}\selectfont TMD}}$" + r", $\gamma + p \rightarrow c + p$" + "\n"
            + rf"$z = 0.5$, $p_\perp = {p_label}$ GeV",
            transform=ax.transAxes, ha="left", va="top", fontsize=23, linespacing=2.0)
    plt.legend(handles=handles, loc="lower right", fontsize=23, frameon=False)

    plt.tight_layout()
    outname = "../plots/charm_fixed_qp_xpom.pdf"
    plt.savefig(outname, bbox_inches="tight")
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
