import glob
import os
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

# Fragmentation functions D(z, Q^2) for c -> D0 at three scales: BCFY, Kniehl-Kramer and HymnD (68% band of
# the replicas). Reads the grids the same way as src/hymnd_grid.cpp. HymnD replicas: HYMND_DIR or the LHAPDF folder.

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 26,
    "xtick.labelsize": 22.5,
    "ytick.labelsize": 22.5,
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
    "xtick.color": "gray",
    "ytick.color": "gray",
    "xtick.labelcolor": "black",
    "ytick.labelcolor": "black",
    "xtick.major.pad": 8,
    "ytick.major.pad": 8,
})

MC = 1.5
CHARM = 4
BCFY_FILE = "../input/BCFY_EKO/bcfy_eko_0000.dat"
KK_FILE = "../input/KK_EKO/kk_eko_0000.dat"
HYMND_SET = "prompt-D0-1-109"
PANELS = [
    (MC, r"$Q_0=m_c=1.5$ GeV", None),
    (np.sqrt(3.0**2 + MC**2), r"$Q=m_t={:.2f}$ GeV".format(np.sqrt(3.0**2 + MC**2)), r"$p_{D\perp}=3$ GeV"),
    (np.sqrt(9.0**2 + MC**2), r"$Q=m_t={:.2f}$ GeV".format(np.sqrt(9.0**2 + MC**2)), r"$p_{D\perp}=9$ GeV"),
]
COLORS = {"BCFY": "#2166ac", "KK": "#b2182b", "HymnD": "#1a7f4b"}


def hymnd_dir():
    """Folder with all HymnD members: HYMND_DIR, else the LHAPDF data folder, else ../input/HymnD."""
    if os.environ.get("HYMND_DIR"):
        return os.environ["HYMND_DIR"]
    try:
        datadir = subprocess.run(["lhapdf-config", "--datadir"], capture_output=True, text=True).stdout.strip()
    except OSError:
        datadir = ""
    for d in [os.path.join(datadir, HYMND_SET), os.path.expanduser(f"~/.local/share/LHAPDF/{HYMND_SET}")]:
        if datadir and len(glob.glob(f"{d}/{HYMND_SET}_*.dat")) > 1:
            return d
    return f"../input/HymnD"


def read_blocks(filename, flavor=CHARM):
    """Q blocks of an LHAPDF grid file: list of (x, Q, values[ix, iq]) for one flavor."""
    with open(filename) as f:
        text = f.read()
    blocks = []
    for part in text.split("---")[1:]:
        lines = part.strip().split("\n")
        if len(lines) < 4:
            continue
        x = np.array(lines[0].split(), dtype=float)
        q = np.array(lines[1].split(), dtype=float)
        col = [int(v) for v in lines[2].split()].index(flavor)
        rows = np.array(" ".join(lines[3:3 + len(x) * len(q)]).split(), dtype=float)
        values = rows.reshape(len(x) * len(q), -1)[:, col].reshape(len(x), len(q))
        blocks.append((x, q, values))
    return blocks


def D_at(blocks, Q):
    """D(z) at scale Q: linear interpolation in Q, as in src/hymnd_grid.cpp. Returns (z, D)."""
    Q = min(max(Q, blocks[0][1][0]), blocks[-1][1][-1])
    for x, q, values in blocks:
        if q[0] <= Q <= q[-1]:
            hi = max(1, min(np.searchsorted(q, Q), len(q) - 1))
            t = (Q - q[hi - 1]) / (q[hi] - q[hi - 1])
            xf = values[:, hi - 1] * (1 - t) + values[:, hi] * t
            return x, xf / x
    raise ValueError(f"Q={Q} not in the grid")


def main():
    bcfy = read_blocks(BCFY_FILE)
    kk = read_blocks(KK_FILE)
    hdir = hymnd_dir()
    members = sorted(glob.glob(f"{hdir}/{HYMND_SET}_*.dat"))
    if not members:
        raise SystemExit(f"No HymnD files found in {hdir}.")
    hymnd = [read_blocks(m) for m in members]
    band = len(hymnd) > 2
    print(f"HymnD: {len(hymnd)} members from {hdir}")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.6), sharey=True, gridspec_kw={"wspace": 0})
    for i, (ax, (Q, label, sublabel)) in enumerate(zip(axes, PANELS)):
        z, d = D_at(bcfy, Q)
        ax.plot(z, d, color=COLORS["BCFY"], lw=2)
        z, d = D_at(kk, Q)
        ax.plot(z, d, color=COLORS["KK"], lw=2, ls="--")
        z, central = D_at(hymnd[0], Q)
        if band:
            replicas = np.array([D_at(m, Q)[1] for m in hymnd[1:]])
            lo, hi = np.percentile(replicas, [16, 84], axis=0)
            ax.fill_between(z, lo, hi, color=COLORS["HymnD"], alpha=0.25, lw=0)
        ax.plot(z, central, color=COLORS["HymnD"], lw=2, ls="-.")
        text = label if sublabel is None else label + "\n" + sublabel
        ax.text(0.08, 0.90, text, transform=ax.transAxes, ha="left", va="top", fontsize=21, linespacing=1.65)
        ax.set_xlim(0, 1)
        ax.set_xlabel(r"$z_h$", labelpad=6)
        ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        hidden = ([] if i == 0 else [0]) + ([] if i == len(axes) - 1 else [1])
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _, h=hidden: "" if round(v, 6) in h else f"{v:g}"))

    axes[0].set_ylim(0, 7.5)
    axes[0].set_ylabel(r"$D_{c\to D^0}(z_h, Q^2)$", labelpad=16)

    handles = [
        Line2D([0], [0], color=COLORS["BCFY"], lw=3, label="BCFY"),
        Line2D([0], [0], color=COLORS["KK"], lw=3, ls="--", label=r"Kniehl \& Kramer"),
        (Patch(facecolor=COLORS["HymnD"], alpha=0.25, lw=0), Line2D([0], [0], color=COLORS["HymnD"], lw=3, ls="-."))
        if band else Line2D([0], [0], color=COLORS["HymnD"], lw=3, ls="-."),
    ]
    labels = ["BCFY", r"Kniehl \& Kramer", r"HymnD (68\% CL)" if band else "HymnD"]
    axes[0].legend(handles, labels, loc="center left", bbox_to_anchor=(0.02, 0.62), fontsize=20, frameon=False)

    outname = "../plots/fragmentation_comparison.pdf"
    plt.savefig(outname, bbox_inches="tight")
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
