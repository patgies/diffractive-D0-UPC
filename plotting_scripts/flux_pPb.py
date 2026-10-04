import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

from flux_comparison import load_scan, CHANNEL_COLORS   # also sets the plot style

# Photon flux of the lead nucleus in p+Pb at 8.16 TeV (as Fig. 11 of arXiv:2606.05469):
# "AnAn" is the Woods-Saxon flux with Gamma_pA(b), "PL(AnAn)" the point-like flux with b_min = 1.1 * 7.1 fm.
# Inputs: output/flux_scan/flux_scan_pPb_{WS,PL}.dat, made by run_scripts/run_flux_scan.sh.

SCAN_DIR = "../output/flux_scan"
CURVES = [   # AnAn first, so the dotted PL curve is drawn on top of it
    ("WS", "AnAn", r"$AnAn$", "-"),
    ("PL", "PL(AnAn)", r"PL$(AnAn)$", ":"),
]


def main():
    scans = {model: load_scan(f"{SCAN_DIR}/flux_scan_pPb_{model}.dat")["pA"] for model, _, _, _ in CURVES}

    fig, ax = plt.subplots(figsize=(8, 7))
    for model, channel, label, linestyle in CURVES:
        z, f = scans[model]
        ax.plot(z, z * f, color=CHANNEL_COLORS[channel], linestyle=linestyle, lw=3, label=label)

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylim(1e-3, 4e2); ax.set_xlim(1e-4, 0.3)
    ax.set_xlabel(r"$z_\gamma$", labelpad=6)
    ax.set_ylabel(r"$z_\gamma\, dF/dz_\gamma$", labelpad=8)
    ax.text(0.95, 0.93, "p-Pb 8.16 TeV", transform=ax.transAxes, ha="right", va="top", fontsize=26)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(LogLocator(base=10.0, numticks=20))
        axis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20))

    ax.legend(loc="lower left", bbox_to_anchor=(0.02, 0.03), fontsize=24, labelspacing=0.9, frameon=False)
    ax.tick_params(labelsize=26)
    ax.tick_params(axis="x", pad=9)
    ax.tick_params(axis="y", pad=6)

    plt.tight_layout()
    out = "../plots/flux_pPb.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    print("Saved:", out)


if __name__ == "__main__":
    main()
