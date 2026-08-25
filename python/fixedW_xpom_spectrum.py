import glob
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

# Plots the fixed-W (no photon flux, proton target, no fragmentation)
# diffractive cross section as a function of ln(x_po), for each rapidity
# and K value, from ../files/D0_fixedW_xpom_y*_K*.dat (produced by
# ../run_fixedW_xpom.sh). Same style/layout as xpom_spectrum.py, but for
# the fixed-W calculation (see notes there on why fixedW has no
# LHAPDF/BCFY split or nucleus target).

plt.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{amssymb}",  # for \mathbb{P} below
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


def read_header(filename):
    """Returns (y, K) read from the '# fixed rapidity y :' / '# fixed K :' header lines."""
    y = None
    K = None
    with open(filename) as f:
        for line in f:
            if "fixed rapidity y" in line:
                y = float(line.split(":")[-1])
            elif "fixed K" in line:
                K = float(line.split(":")[-1])
    if y is None or K is None:
        raise ValueError(f"missing rapidity/K header in {filename}")
    return y, K


def load_results():
    """Read all files and return {(y, K): [(x_po, dsigma_dlnxpo), ...]}."""
    pattern = "../files/D0_fixedW_xpom_y*_K*.dat"
    results = {}
    for filename in sorted(glob.glob(pattern)):
        y, K = read_header(filename)
        data = np.loadtxt(filename, comments='#')
        xpo = data[:, 0]
        diff = data[:, 2]
        # dln(x_po) = dx_po / x_po
        dsigma_dlnxpo = xpo * diff
        results[(y, K)] = sorted(zip(xpo, dsigma_dlnxpo))
    return results


def main():
    results = load_results()
    if not results:
        sys.exit("No files found matching ../files/D0_fixedW_xpom_y*_K*.dat -- "
                  "run ../run_fixedW_xpom.sh first.")

    # K values to include in the plot:
    K_TO_PLOT = [2.0]
    filtered_results = {(y, K): pts for (y, K), pts in results.items() if K in K_TO_PLOT}
    results = filtered_results
    if not results:
        sys.exit(f"No data for K in {K_TO_PLOT} -- check ../files/D0_fixedW_xpom_y*_K*.dat.")

    y_values = sorted(set(y for y, K in results))
    K_values = sorted(set(K for y, K in results))

    blue_ramp = ["#cde2fb", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
    colors = {}
    for i, y in enumerate(y_values):
        step = round(i * (len(blue_ramp) - 1) / max(len(y_values) - 1, 1))
        colors[y] = blue_ramp[step]

    linestyle_cycle = ['-', '--', ':', '-.']
    y_linestyles = {}
    for i, y in enumerate(y_values):
        y_linestyles[y] = linestyle_cycle[i % len(linestyle_cycle)]

    plt.figure(figsize=(7.5, 6.5))
    for (y, K), points in sorted(results.items()):
        xpo_values = [p[0] for p in points if p[1] > 0]
        dsigma_values = [p[1] for p in points if p[1] > 0]
        if not xpo_values:
            continue
        ln_xpo = np.log(xpo_values)
        plt.plot(ln_xpo, dsigma_values, color=colors[y], linestyle=y_linestyles[y])

    plt.yscale("log")
    plt.xlabel(r"$\ln(x_{\mathbb{P}})$")
    plt.ylabel(r"$d\sigma/dp_{c\perp}dy\,d\ln(x_{\mathbb{P}})$ [mb/GeV]", labelpad=15)

    # Secondary top axis showing x_po
    ax = plt.gca()
    ax.tick_params(axis='x', which='both', top=False)
    secax = ax.secondary_xaxis('top', functions=(np.exp, np.log))
    secax.set_xscale('log')
    secax.set_xlabel(r"$x_{\mathbb{P}}$", labelpad=10)

    ax.set_title(r"Diffractive charm with proton target, fixed $q^+=3p^+$", pad=15)

    K_labels = [f"$p_{{c\\perp}}$={K:g} GeV" for K in K_values]
    K_title = ", ".join(K_labels)
    y_handles = [Line2D([0], [0], color=colors[y], linestyle=y_linestyles[y], label=f"y={y:g}") for y in reversed(y_values)]
    plt.legend(handles=y_handles, loc="upper left", fontsize=13, title=K_title, title_fontsize=13)

    plt.tight_layout()
    outname = "../plots/fixedW_xpom_spectrum.pdf"
    plt.savefig(outname, dpi=150)
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
