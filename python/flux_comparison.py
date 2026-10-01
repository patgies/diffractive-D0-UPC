import csv
import os
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, FixedLocator, FuncFormatter, NullFormatter

# make the plot look nicer (same style as cross_section.py / xpom_spectrum.py / fixed_qp_spectrum.py)
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
# The effective fluxes use sigma_NN = 90.8533 mb (inputs/Gamma_AA_sigma90.85.dat),
# the value behind the Starlight table (see its header: sqrt(s_NN) = 5.36 TeV);
# the physics pipeline uses 92 mb (inputs/Gamma_AA.dat).
# Regenerate with, from the repository root:
#   ./build/bin/scan_flux > out/flux_scan.csv
#   for ch in AnAn An0n Xn0n; do
#     GAMMA_AA_FILE=./inputs/Gamma_AA_sigma90.85.dat ./build/bin/scan_flux_eff $ch > out/flux_scan_eff_$ch.csv
#   done

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
# identity; line style (TABLE_STYLE/OURS_STYLE below) carries the model.
PAIRS = [
    ("PL(AnAn)", PL, "AnAn", r"PL$(AnAn)$"),
    ("AnAn",     WS, "AnAn", r"$AnAn$"),
    ("An0n",     WS, "An0n", r"$An0n$"),
    XN0N_ENTRY,
]
# dataviz skill categorical slots 1-4, fixed order (adjacent-pair CVD deltaE
# 9.1/8.4 light/dark, normal-vision 19.6/19.3 -- validated for the line-chart
# adjacent pairlist, not the stricter all-pairs one).
# PL(AnAn) is the point-like baseline, not a nuclear channel -> neutral ink,
# so it doesn't compete with (and get confused with) the AnAn hue it overlaps.
CHANNEL_COLORS = {
    "PL(AnAn)": "#333333",  # neutral baseline
    "AnAn":     "#d62728",  # red
    "An0n":     "#2563b8",  # medium-dark blue
    "Xn0n":     "#1a7f4b",  # dark green
}
# Starlight (table) = thick pale band underneath; P_b (this code) = thin
# short-dashed line on top. The two agree almost everywhere, so a same-width
# solid/dashed pair just hides the dashes under the solid line.
TABLE_STYLE = dict(ls="-", lw=3.5, alpha=0.35, solid_capstyle="butt")
OURS_STYLE  = dict(ls=(0, (3, 1.5)), lw=1.5)

if __name__ == "__main__":
    TABLED = [p for p in PAIRS if p[1] is not None]   # Xn0n has no table to compare
    print(f"{'y':>9} " + " ".join(f"{c[3][:14]:>16}" for c in TABLED) + "   (ours/table)")
    for yy in [1e-4, 1e-3, 1e-2, 1e-1, 0.5]:
        line = f"{yy:9.1e} "
        for ch, tab, key, _label in TABLED:
            y, f = scan[ch]
            # Interpolate the ratio (smooth, ~1), not the flux itself: the
            # scan grid is too coarse for log-log interpolation of the
            # steeply falling flux at large y, which fakes a deficit
            # (e.g. 0.84 instead of 1.03 for PL(AnAn) at y=0.5).
            ratio = f / interp_flux(tab, key, y)
            line += f"{np.interp(np.log(yy), np.log(y), ratio):16.3f} "
        print(line)

    YMIN = 1e-3   # lower edge of the upper panel
    # Ratio lines in the same solid channel colors as the upper panel.
    RATIO_STYLE = {"PL(AnAn)": dict(lw=1.5),
                   "AnAn":     dict(lw=1.5),
                   "An0n":     dict(lw=1.5)}
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(7.5, 7.4), sharex=True, gridspec_kw={"height_ratios": [3, 0.75]})
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
            ax.plot(y, y * f, color=color, **OURS_STYLE)
            continue
        # Lower panel: P_b / Starlight, i.e. how well this code reproduces
        # the tabulated flux, for every channel (including PL(AnAn)).
        # Evaluated only at the table's own Chebyshev nodes, where the
        # reference is exact: in between, its barycentric interpolation
        # (Eq. 34 of arXiv:2404.09731) is off by up to ~0.5% at z > 0.05,
        # which showed up as a spurious wiggle. Our flux is densely sampled
        # and smooth, so a cubic spline in (log z, log f) is exact enough.
        ok = (f > 0) & (y < 0.7)   # PL scan underflows to 0 at y > 0.8
        ours = CubicSpline(np.log(y[ok]), np.log(f[ok]))
        nodes = tab["y4"] ** 4
        at = (nodes >= y[ok][0]) & (nodes <= y[ok][-1])
        axr.plot(nodes[at], np.exp(ours(np.log(nodes[at])) - tab[key][at]),
                 color=color, ls="-", **RATIO_STYLE[ch])
        ax.plot(yy, yy * interp_flux(tab, key, yy), color=color, **TABLE_STYLE)
        ax.plot(y, y * f, color=color, **OURS_STYLE)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylim(YMIN, 4e2); ax.set_xlim(1e-4, 0.12)
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
    channel_handles = [Line2D([], [], color=CHANNEL_COLORS[ch], lw=3, ls="-", label=label) for ch, _, _, label in PAIRS]
    source_handles = [Line2D([], [], color="0.3", **OURS_STYLE, label=r"$P_b$ (Baur \textit{et al.})"),
                       Line2D([], [], color="0.3", **TABLE_STYLE, label=r"Starlight (Eskola \textit{et al.})")]
    channel_legend = ax.legend(handles=channel_handles, loc="lower left", bbox_to_anchor=(0.02, 0.0),
                                fontsize=15, title="Channel", title_fontsize=15, frameon=False)
    ax.add_artist(channel_legend)
    ax.legend(handles=source_handles, loc="lower left", bbox_to_anchor=(0.3, 0.0),
              fontsize=15, title="Model", title_fontsize=15, frameon=False)
    # Ratio-panel legend: one entry per ratio line, drawn in its ratio style.
    ratio_handles = [Line2D([], [], color=CHANNEL_COLORS[ch], ls="-", **RATIO_STYLE[ch], label=label)
                     for ch, tab, _, label in PAIRS if tab is not None]
    axr.legend(handles=ratio_handles, loc="lower left", ncol=3, fontsize=15,
               frameon=False, handlelength=1.8, columnspacing=1.2)
    axr.set_xscale("log"); axr.set_ylabel(r"$P_b$/Starlight", labelpad=15, fontsize=17); axr.set_xlabel(r"$z_\gamma$", labelpad=6)
    # x stops at z = 0.12 (beyond, flux < 1e-3 and the ratio is numerical
    # noise), so the ratio stays below ~1.25; y from 0.95 to 1.05 (AnAn/An0n run
    # off the top at z ~ 0.06): linear scale.
    axr.set_yscale("linear"); axr.set_ylim(0.95, 1.05)
    axr.yaxis.set_major_locator(FixedLocator([0.95, 1.0, 1.05]))
    axr.yaxis.set_minor_locator(FixedLocator(np.arange(0.95, 1.051, 0.01)))
    axr.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    axr.yaxis.set_minor_formatter(NullFormatter())
    fig.align_ylabels([ax, axr])
    plt.tight_layout()
    out = "../plots/flux_comparison.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    print("Saved:", out)
