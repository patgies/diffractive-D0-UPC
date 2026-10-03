import glob
import math
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.integrate import simpson
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))   # alphas_running.py
from alphas_running import alphas_run


plt.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{amssymb}",
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

# Reads the D0_xpom output of run_xpom.sh, integrates over b_d and plots
# the diffractive cross section versus ln(x_po).

alphae = 1/137
mc     = 1.5
e_c    = 2/3
Nc     = 3

FMGEV = 5.068
GEVSQR_TO_NB = 1.0e7 / (FMGEV * FMGEV)
GEVSQR_TO_MB = GEVSQR_TO_NB * 1e-6

NUCLEUS = os.environ.get("NUCLEUS", "Pb")
FRAG    = os.environ.get("FRAG_TYPE", "BCFY")
CHANNEL = os.environ.get("CHANNEL", "An0n").translate(str.maketrans('', '', '() '))


def read_header(filename):
    """Returns (y, pD0) read from the '# fixed rapidity y :' / '# fixed pD0 :' header lines."""
    y = None
    pt = None
    with open(filename) as f:
        for line in f:
            if "fixed rapidity y" in line:
                y = float(line.split(":")[-1])
            elif "fixed pD0" in line:
                pt = float(line.split(":")[-1])
    if y is None or pt is None:
        raise ValueError(f"missing rapidity/pD0 header in {filename}")
    return y, pt


def read_data_file(filename):
    """Returns three plain lists (b_d, x_po, dsigma), skipping comment lines."""
    b_d_list = []
    xpo_list = []
    dsigma_list = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            b_d, xpo, dsigma = line.split()
            b_d_list.append(float(b_d))
            xpo_list.append(float(xpo))
            dsigma_list.append(float(dsigma))
    return b_d_list, xpo_list, dsigma_list


def group_by_xpo(b_d_list, xpo_list, dsigma_list):
    groups = {}
    for b_d, xpo, dsigma in zip(b_d_list, xpo_list, dsigma_list):
        if xpo not in groups:
            groups[xpo] = []
        groups[xpo].append((b_d, dsigma))
    return groups


def integrate_over_b_d(pairs):
    """Integral of b_d*dsigma(b_d) db_d with Simpson's rule."""
    pairs = sorted(pairs)
    if pairs[0][0] < 0:
        raise ValueError(
            "Negative b_d found -- datafile wasn't a Glauber sample "
            "(e.g. the plain proton dipole has no b_d to integrate over)."
        )
    b_d_values = [pair[0] for pair in pairs]
    weighted = [b_d * dsigma for b_d, dsigma in pairs]
    return simpson(weighted, x=b_d_values)


def diffractive_prefactor(pt):
    """Prefactor of the cross section. No sigma0: the b_d integral already gives the area."""
    alphas = alphas_run(math.sqrt(pt**2 + mc**2))
    return alphas * alphae * e_c**2 * (Nc**2 - 1) / (8 * math.pi**4)


def exclusive_prefactor():
    """Prefactor for exclusive. It does not depend on pt or alpha_s. No sigma0."""
    return alphae * Nc * e_c**2 / (2 * math.pi**2)


def load_results(process="diffractive"):
    """Read all files for one process ("exclusive"/"diffractive") and return
    {(y, pt): [(x_po, dsigma_dy_dpt_dlnxpo), ...]}."""
    prefactor = diffractive_prefactor if process == "diffractive" else (lambda pt: exclusive_prefactor())
    name = f"D0_{process}_xpom_{FRAG}_{CHANNEL}_{NUCLEUS}_y*_pt*.dat"
    xpom_dir = f"../output/{CHANNEL}/xpom"
    results = {}
    for filename in sorted(glob.glob(f"{xpom_dir}/{FRAG}/{name}") or glob.glob(f"{xpom_dir}/{name}")):
        y, pt = read_header(filename)
        b_d_list, xpo_list, dsigma_list = read_data_file(filename)
        xpo_groups = group_by_xpo(b_d_list, xpo_list, dsigma_list)

        points = []
        for xpo, pairs in sorted(xpo_groups.items()):
            b_d_integral = integrate_over_b_d(pairs)
            dsigma_dxpo = (2*math.pi) * b_d_integral * prefactor(pt) \
                          * (2*math.pi) * pt * GEVSQR_TO_MB
            dsigma_dlnxpo = xpo * dsigma_dxpo
            points.append((xpo, dsigma_dlnxpo))
        results[(y, pt)] = sorted(points)
    return results


def main():
    results = load_results("diffractive")
    if not results:
        sys.exit(f"No files found for NUCLEUS={NUCLEUS}, FRAG_TYPE={FRAG}, CHANNEL={CHANNEL} "
                  f"in ../output/{CHANNEL}/xpom/ -- run ../run_scripts/run_xpom.sh first.")

    PT_TO_PLOT = [
        # 1.0,
        2.0,
        # 4.0,
        # 8.0,
    ]
    filtered_results = {}
    for (y, pt), points in results.items():
        if pt in PT_TO_PLOT:
            filtered_results[(y, pt)] = points
    results = filtered_results

    y_values = sorted(set(y for y, pt in results))
    pt_values = sorted(set(pt for y, pt in results))

    palette = ["#2166ac", "#4393c3", "#92c5de", "#f4a582", "#d6604d", "#b2182b"]
    if len(y_values) <= 4:
        palette = ["#2166ac", "#67a9cf", "#ef8a62", "#b2182b"]
    colors = {y: palette[round(i * (len(palette) - 1) / max(len(y_values) - 1, 1))]
              for i, y in enumerate(y_values)}
    linestyle_cycle = ["-", "--", ":", "-.", (0, (5, 1.5, 1, 1.5, 1, 1.5)), (0, (8, 2))]
    y_linestyles = {y: linestyle_cycle[i % len(linestyle_cycle)] for i, y in enumerate(y_values)}

    plt.figure(figsize=(8, 7))
    for (y, pt), points in sorted(results.items()):
        xpo_values = [p[0] for p in points]
        dsigma_values = [p[1] for p in points]
        ln_xpo = np.log(xpo_values)
        plt.plot(ln_xpo, dsigma_values, color=colors[y], linestyle=y_linestyles[y], lw=2.0)

    plt.yscale("log")
    plt.ylim(1e-9, 1e5)
    plt.yticks(10.0**np.arange(-9, 6, 2))
    plt.xlim(-10, -2)
    plt.xlabel(r"$\ln(x_{\mathbb{P}})$", labelpad=14)
    plt.ylabel(r"$d\sigma/dy\,dp_{D^0\perp}\,d\ln(x_{\mathbb{P}})$ [mb/GeV]", labelpad=16)

    ax = plt.gca()
    ax.tick_params(axis='x', which='both', top=False)
    secax = ax.secondary_xaxis('top', functions=(np.exp, np.log))
    secax.set_xscale('log')
    secax.set_xlabel(r"$x_{\mathbb{P}}$", labelpad=23)

    process_label = r"Diffractive$_{\mbox{\fontsize{13}{13}\selectfont TMD}}$"
    pt_label = ", ".join(f"{pt:g}" for pt in pt_values)
    ax.text(0.95, 0.95,
            f"{NUCLEUS}-{NUCLEUS} 5.36 TeV\n{process_label}\n{CHANNEL}, {FRAG}\n"
            f"$p_{{D^0\\perp}} = {pt_label}$ GeV",
            transform=ax.transAxes, ha="right", va="top", fontsize=21, linespacing=1.8)

    y_handles = [Line2D([0], [0], color=colors[y], linestyle=y_linestyles[y], lw=3, label=f"$y={y:g}$")
                 for y in y_values]
    plt.legend(handles=y_handles, loc="upper left", ncol=1, fontsize=21, frameon=False)

    plt.tight_layout()
    outname = f"../plots/xpom_{CHANNEL}_{NUCLEUS}{NUCLEUS}.pdf"
    plt.savefig(outname, dpi=150)
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
