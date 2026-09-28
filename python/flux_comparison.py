import csv
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, FixedLocator, NullFormatter, FuncFormatter

# make the plot look nicer (same style as cross_section.py / xpom_spectrum.py / fixedW_spectrum.py)
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
    "xtick.color": "0.4",
    "ytick.color": "0.4",
    "xtick.labelcolor": "black",
    "ytick.labelcolor": "black",
    "xtick.major.pad": 8,
    "ytick.major.pad": 4,
})

# Compares the photon flux of this code with the tabulated effective fluxes
# of arXiv:2404.09731 (inputs/photon_flux/log-flux-tbl-*.dta, P. Paakkinen),
# evaluated with the paper's barycentric Lagrange interpolation on Chebyshev
# nodes in y^(1/4) (Eqs. 34-37). AnAn/An0n/Xn0n use this code's *effective*
# flux (src/scan_flux_eff.cpp -> out/flux_scan_eff_<channel>.csv, Eq. 4 of
# arXiv:2404.09731 -- see photon_flux_discrepancy.pdf for why the simple
# single-b treatment used to live here instead, and why it undershoots the
# table by up to 33x at large y). PL(AnAn) is the deliberate point-like
# baseline, not the effective flux, so it still comes from the simple scan
# (src/scan_flux.cpp -> out/flux_scan.csv).
# Regenerate with, from python/:
#   ../build/bin/scan_flux > ../out/flux_scan.csv
#   ../build/bin/scan_flux_eff AnAn > ../out/flux_scan_eff_AnAn.csv
#   ../build/bin/scan_flux_eff An0n > ../out/flux_scan_eff_An0n.csv
#   ../build/bin/scan_flux_eff Xn0n > ../out/flux_scan_eff_Xn0n.csv

TABLE_DIR = "../inputs/photon_flux"
OUT_DIR = "../out"


def load_table(path):
    rows = [l.split() for l in open(path) if l.strip() and not l.startswith("#")]
    a = np.array(rows, dtype=float)
    return {"y4": a[:, 2], "AnAn": a[:, 3], "An0n": a[:, 4]}


def interp_flux(tab, key, y):
    """f(y) = exp[ sum_j w_j/(y4-y4_j) log f_j / sum_j w_j/(y4-y4_j) ]  (Eq. 34)."""
    y4j, logf = tab["y4"], tab[key]
    n = len(y4j) - 1
    w = (-1.0) ** np.arange(n + 1)
    w[0] *= 0.5
    w[-1] *= 0.5
    out = []
    for yy in np.atleast_1d(y):
        d = yy ** 0.25 - y4j
        k = np.argmin(abs(d))
        if abs(d[k]) < 1e-14:
            out.append(np.exp(logf[k]))
            continue
        r = w / d
        out.append(np.exp(np.sum(r * logf) / np.sum(r)))
    return np.array(out)


def load_csv(fn):
    d = {}
    for r in csv.DictReader(open(fn)):
        d.setdefault(r["channel"], []).append((float(r["y"]), float(r["dN_dy"])))
    return {k: (np.array([p[0] for p in v]), np.array([p[1] for p in v])) for k, v in d.items()}


PL = load_table(os.path.join(TABLE_DIR, "log-flux-tbl-PL.dta"))
WS = load_table(os.path.join(TABLE_DIR, "log-flux-tbl-WS.dta"))
scan = {"PL(AnAn)": load_csv(os.path.join(OUT_DIR, "flux_scan.csv"))["PL(AnAn)"]}
for ch in ["AnAn", "An0n", "Xn0n"]:
    scan[ch] = load_csv(os.path.join(OUT_DIR, f"flux_scan_eff_{ch}.csv"))[ch]

# Xn0n has no column in Petja's tables at all (only AnAn/An0n) -- plotted
# as a mine-only curve (no table comparison, no ratio-panel line).
XN0N_ENTRY = ("Xn0n", None, None, r"$Xn0n$")

# (our channel, table, table key, label) -- color now carries the channel
# identity; linestyle (TABLE_LS/OURS_LS below) carries the model.
PAIRS = [
    ("PL(AnAn)", PL, "AnAn", r"PL$(AnAn)$"),
    ("AnAn",     WS, "AnAn", r"$AnAn$"),
    ("An0n",     WS, "An0n", r"$An0n$"),
    XN0N_ENTRY,
]
# dataviz skill categorical slots 1-4, fixed order (adjacent-pair CVD deltaE
# 9.1/8.4 light/dark, normal-vision 19.6/19.3 -- validated for the line-chart
# adjacent pairlist, not the stricter all-pairs one).
CHANNEL_COLORS = {
    "PL(AnAn)": "#2a78d6",  # slot 1, blue
    "AnAn":     "#eb6834",  # slot 2, orange
    "An0n":     "#1baf7a",  # slot 3, aqua
    "Xn0n":     "#eda100",  # slot 4, yellow
}
TABLE_LS = "-"   # Starlight (table)
OURS_LS  = "--"  # P_b (this code)

if __name__ == "__main__":
    TABLED = [p for p in PAIRS if p[1] is not None]   # Xn0n has no table to compare
    print(f"{'y':>9} " + " ".join(f"{c[3][:14]:>16}" for c in TABLED) + "   (ours/table)")
    for yy in [1e-4, 1e-3, 1e-2, 1e-1, 0.5]:
        line = f"{yy:9.1e} "
        for ch, tab, key, _label in TABLED:
            y, f = scan[ch]
            ours = np.interp(np.log(yy), np.log(y), np.log(f))
            line += f"{np.exp(ours) / interp_flux(tab, key, yy)[0]:16.3f} "
        print(line)

    fig, (ax, axr) = plt.subplots(2, 1, figsize=(7.5, 7), sharex=True, gridspec_kw={"height_ratios": [6, 1]})
    yy = np.logspace(-4, 0, 400)
    # bmax = 60/(y*m_N) (as in the old DiffractiveD0/testFlux code, and now
    # also what main.cpp/main_xpom.cpp use, worst-case-adjusted for the
    # (pD0, y) they're run at -- see src/scan_flux.cpp and main.cpp) --
    # large enough at every y to fully capture the b->infinity flux integral.
    for ch, tab, key, label in PAIRS:
        y, f = scan[ch]
        color = CHANNEL_COLORS[ch]
        if tab is None:
            # No table for Xn0n at all: mine-only, no ratio-panel line.
            ax.plot(y, y * f, color=color, ls=OURS_LS, lw=1.5)
            continue
        # Lower panel: P_b / Starlight, i.e. how well this code reproduces
        # the tabulated flux, for every channel (including PL(AnAn)).
        axr.plot(y, f / interp_flux(tab, key, y), color=color, ls="-", lw=1.5)
        ax.plot(yy, yy * interp_flux(tab, key, yy), color=color, ls=TABLE_LS, lw=2)
        ax.plot(y, y * f, color=color, ls=OURS_LS, lw=1.5)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylim(1e-3, 4e2); ax.set_xlim(1e-4, 0.5)
    ax.set_ylabel(r"$z_\gamma\, dF/dz_\gamma$", labelpad=8)
    ax.text(0.95, 0.95, "Pb-Pb 5.36 TeV", transform=ax.transAxes, ha="right", va="top", fontsize=16)
    for a in (ax, axr):
        a.xaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
        a.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20))
        a.yaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
        a.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20))

    # Upper axis: photon energy omega = y*sqrt(s_NN)/2 [GeV] (pure rescaling
    # of y, same sqrt(s_NN)=5360 GeV used in src/scan_flux.cpp).
    SQRT_S = 5360.0
    y_to_omega = lambda y: y * SQRT_S / 2.0
    omega_to_y = lambda w: w * 2.0 / SQRT_S
    secax = ax.secondary_xaxis("top", functions=(y_to_omega, omega_to_y))
    secax.set_xlabel(r"$\omega$ [GeV]", labelpad=16)
    secax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    secax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20))
    secax.tick_params(axis="x", which="major", pad=4)

    # Two separate legends, stacked in the lower-left corner (channel above
    # source -- the citation labels are too long to sit side by side at
    # this font size).
    channel_handles = [Line2D([], [], color=CHANNEL_COLORS[ch], lw=2, ls="-", label=label) for ch, _, _, label in PAIRS]
    source_handles = [Line2D([], [], color="k", lw=2, ls=TABLE_LS, label=r"Starlight (Eskola \textit{et al.})"),
                       Line2D([], [], color="k", lw=2, ls=OURS_LS, label=r"$P_b$ (Baur \textit{et al.})")]
    channel_legend = ax.legend(handles=channel_handles, loc="lower left", bbox_to_anchor=(0.02, 0.22),
                                fontsize=13, title="Channel", title_fontsize=13, frameon=False)
    ax.add_artist(channel_legend)
    ax.legend(handles=source_handles, loc="lower left", bbox_to_anchor=(0.02, 0.0),
              fontsize=13, title="Model", title_fontsize=13, frameon=False)
    axr.set_xscale("log"); axr.set_ylabel(r"$P_b$/Starlight", labelpad=15, fontsize=15); axr.set_xlabel(r"$z_\gamma$", labelpad=6)
    axr.set_yscale("log"); axr.set_ylim(5e-1, 3.0)
    # Plain-decimal ticks on the ratio panel: the default log minor-tick
    # labels (2x,3x,4x,6x...) crowd and overlap over a <1-decade range.
    axr.yaxis.set_major_locator(FixedLocator([0.5, 1, 2]))
    axr.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    axr.yaxis.set_minor_formatter(NullFormatter())
    fig.align_ylabels([ax, axr])
    plt.tight_layout()
    out = "../plots/flux_comparison.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    print("Saved:", out)
