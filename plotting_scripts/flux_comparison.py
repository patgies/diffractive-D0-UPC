import os
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, FixedLocator, FuncFormatter, NullFormatter

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
    "xtick.color": "black",
    "ytick.color": "black",
    "xtick.labelcolor": "black",
    "ytick.labelcolor": "black",
    "xtick.major.pad": 8,
    "ytick.major.pad": 4,
})

# Compares the photon flux of this code with the tables of arXiv:2404.09731. Inputs:
# output/flux_scan/, made by run_scripts/run_flux_scan.sh (sigma_NN = 90.85 mb of the Starlight tables).

TABLE_DIR = "../input/Starlight_photon_flux"
SCAN_DIR = "../output/flux_scan"


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


def load_scan(fn):
    """scan_flux output: "# channel  y  dN_dy" header, one row per (channel, y)."""
    d = {}
    for line in open(fn):
        if not line.strip() or line.startswith("#"):
            continue
        ch, y, dn_dy = line.split()
        d.setdefault(ch, []).append((float(y), float(dn_dy)))
    return {k: (np.array([p[0] for p in v]), np.array([p[1] for p in v])) for k, v in d.items()}


PL = load_table(os.path.join(TABLE_DIR, "log-flux-tbl-PL.dta"))
WS = load_table(os.path.join(TABLE_DIR, "log-flux-tbl-WS.dta"))
scan = {"PL(AnAn)": load_scan(os.path.join(SCAN_DIR, "flux_scan_PL.dat"))["PL(AnAn)"]}
eff = load_scan(os.path.join(SCAN_DIR, "flux_scan_EFF.dat"))
for ch in ["AnAn", "An0n", "Xn0n", "0n0n"]:
    scan[ch] = eff[ch]

XN0N_ENTRY = ("Xn0n", None, None, r"$Xn0n$")
ZN0N_ENTRY = ("0n0n", None, None, r"$0n0n$")

PAIRS = [
    ("PL(AnAn)", PL, "AnAn", r"PL$(AnAn)$"),
    ("AnAn",     WS, "AnAn", r"$AnAn$"),
    ("An0n",     WS, "An0n", r"$An0n$"),
    ZN0N_ENTRY,
    XN0N_ENTRY,
]
CHANNEL_COLORS = {
    "PL(AnAn)": "#333333",
    "AnAn":     "#d62728",
    "An0n":     "#2563b8",
    "Xn0n":     "#1f907e",
    "0n0n":     "#e08a00",
}
TABLE_STYLE = dict(ls="-", lw=6.5, alpha=0.35, solid_capstyle="butt")
OURS_STYLE  = dict(ls=(0, (3, 1.5)), lw=2.0)

if __name__ == "__main__":
    TABLED = [p for p in PAIRS if p[1] is not None]
    print(f"{'y':>9} " + " ".join(f"{c[3][:14]:>16}" for c in TABLED) + "   (ours/table)")
    for yy in [1e-4, 1e-3, 1e-2, 1e-1, 0.5]:
        line = f"{yy:9.1e} "
        for ch, tab, key, _label in TABLED:
            y, f = scan[ch]
            ratio = f / interp_flux(tab, key, y)
            line += f"{np.interp(np.log(yy), np.log(y), ratio):16.3f} "
        print(line)

    YMIN = 1e-3
    RATIO_STYLE = {"AnAn":     dict(lw=1.5),
                   "An0n":     dict(lw=1.5)}
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(8, 9.5), sharex=True, gridspec_kw={"height_ratios": [4, 1]})
    yy = np.logspace(-4, 0, 400)
    for ch, tab, key, label in PAIRS:
        y, f = scan[ch]
        color = CHANNEL_COLORS[ch]
        if tab is None:
            ax.plot(y, y * f, color=color, **OURS_STYLE)
            continue
        ok = (f > 0) & (y < 0.7)
        ours = CubicSpline(np.log(y[ok]), np.log(f[ok]))
        nodes = tab["y4"] ** 4
        at = (nodes >= y[ok][0]) & (nodes <= y[ok][-1])
        if ch in RATIO_STYLE:
            axr.plot(nodes[at], np.exp(ours(np.log(nodes[at])) - tab[key][at]),
                     color=color, ls="-", **RATIO_STYLE[ch])
        ax.plot(yy, yy * interp_flux(tab, key, yy), color=color, **TABLE_STYLE)
        ax.plot(y, y * f, color=color, **OURS_STYLE)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylim(YMIN, 4e2); ax.set_xlim(1e-4, 0.1)
    ax.set_ylabel(r"$z_\gamma\, dF/dz_\gamma$", labelpad=8)
    ax.text(0.95, 0.93, "Pb-Pb 5.36 TeV", transform=ax.transAxes, ha="right", va="top", fontsize=26)
    for a in (ax, axr):
        a.tick_params(labelsize=26)
        a.tick_params(axis="x", pad=9)
        a.tick_params(axis="y", pad=6)
        a.xaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
        a.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20))
        a.yaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
        a.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20))

    SQRT_S = 5360.0
    y_to_omega = lambda y: y * SQRT_S / 2.0
    omega_to_y = lambda w: w * 2.0 / SQRT_S
    secax = ax.secondary_xaxis("top", functions=(y_to_omega, omega_to_y))
    secax.set_xlabel(r"$\omega$ [GeV]", labelpad=26)
    secax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    secax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20))
    secax.tick_params(axis="x", which="major", pad=6, labelsize=26)

    channel_handles = [Line2D([], [], color=CHANNEL_COLORS[ch], lw=3, ls="-", label=label) for ch, _, _, label in PAIRS]
    source_handles = [Line2D([], [], color="0.3", **OURS_STYLE, label=r"$P_b$ (Baur \textit{et al.})"),
                       Line2D([], [], color="0.3", **TABLE_STYLE, label=r"Starlight (Eskola \textit{et al.})")]
    channel_legend = ax.legend(handles=channel_handles, loc="lower left", bbox_to_anchor=(0.02, 0.215),
                                fontsize=18, labelspacing=0.55, frameon=False)
    ax.add_artist(channel_legend)
    ax.legend(handles=source_handles, loc="lower left", bbox_to_anchor=(0.02, 0.0),
              fontsize=18, labelspacing=0.55, title="Model", title_fontsize=18, frameon=False)
    ratio_handles = [Line2D([], [], color=CHANNEL_COLORS[ch], ls="-", **RATIO_STYLE[ch], label=label)
                     for ch, _, _, label in PAIRS if ch in RATIO_STYLE]
    axr.legend(handles=ratio_handles, loc="upper left", bbox_to_anchor=(0.0, 1.04), ncol=3, fontsize=20,
               frameon=False, handlelength=1.8, columnspacing=1.2)
    axr.set_xscale("log"); axr.set_ylabel(r"$P_b$/Starlight", labelpad=15, fontsize=20); axr.set_xlabel(r"$z_\gamma$", labelpad=6)
    axr.set_yscale("linear"); axr.set_ylim(0.98, 1.02)
    axr.yaxis.set_major_locator(FixedLocator([0.98, 1.0, 1.02]))
    axr.yaxis.set_minor_locator(FixedLocator(np.arange(0.98, 1.0201, 0.005)))
    axr.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    axr.yaxis.set_minor_formatter(NullFormatter())
    fig.align_ylabels([ax, axr])
    plt.tight_layout()
    out = "../plots/flux_comparison.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    print("Saved:", out)
