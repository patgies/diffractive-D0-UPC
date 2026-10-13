import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Pb-Pb photon flux z F(z) with our effective flux (EFF, arXiv:2404.09731 Eq. 4: flux at the struck nucleon)
# and with the flux at the centre of the target (FLUX_MODEL=WS, as Eq. (20) of arXiv:2606.05469),
# both with sigma_NN = 92 mb, and the point-like AnAn flux. Lower panel: ratios to F.
# Data: scan_flux EFF/WS AnAn An0n Xn0n with GAMMA_AA_FILE=input/WS_photon_flux/Gamma_AA.dat.
from flux_comparison import load_scan, CHANNEL_COLORS   # also sets the plot style

CHANNELS = ["AnAn", "An0n", "Xn0n"]
EFF_STYLE = dict(ls="-", lw=2.5)
WS_STYLE = dict(ls="--", lw=2.0)
PL_STYLE = dict(ls=":", lw=2.5, color="black")


def main():
    eff = load_scan("../output/flux_scan/flux_scan_EFF_sigma92.dat")
    ws = load_scan("../output/flux_scan/flux_scan_WS_sigma92.dat")
    pl = load_scan("../output/flux_scan/flux_scan_PL.dat")   # point-like, b > b_min = 14.2 fm (no Gamma_AA)

    fig, (ax, ax_ratio) = plt.subplots(2, 1, figsize=(8, 8.1), sharex=True,
                                       gridspec_kw={"height_ratios": [5, 1.1], "hspace": 0.06})
    for ch in CHANNELS:
        z_ws, f_ws = ws[ch]
        for scan, style in [(eff, EFF_STYLE), (ws, WS_STYLE)]:
            z, f = scan[ch]
            ax.plot(z, z * f, color=CHANNEL_COLORS[ch], **style)
        z, f = eff[ch]
        ax_ratio.plot(z, f / np.interp(z, z_ws, f_ws), color=CHANNEL_COLORS[ch], **EFF_STYLE)
    z, f = pl["PL(AnAn)"]
    ax.plot(z, z * f, **PL_STYLE)
    z_ws, f_ws = ws["AnAn"]
    ax_ratio.plot(z, f / np.interp(z, z_ws, f_ws), **PL_STYLE)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-4, 0.1)
    ax.set_ylim(1e-3, 3e2)
    ax.set_ylabel(r"$z_\gamma\,F(z_\gamma)$", labelpad=10)

    ax_ratio.set_yscale("log")
    ax_ratio.set_ylim(0.8, 40)
    ax_ratio.set_yticks([1, 10])
    ax_ratio.set_yticklabels(["1", "10"])
    ax_ratio.set_xlabel(r"$z_\gamma$", labelpad=10)
    ax_ratio.set_ylabel(r"ratio to $F$", labelpad=10, fontsize=20)

    ax.text(0.95, 0.95, r"Pb-Pb 5.36 TeV" + "\n" + r"$\sigma_{NN} = 92$ mb", transform=ax.transAxes,
            ha="right", va="top", fontsize=22, linespacing=1.4)
    channel_handles = [Line2D([0], [0], color=CHANNEL_COLORS[ch], lw=3, label=ch) for ch in CHANNELS]
    # the symbols in a fixed-width box, so that the descriptions start in the same column
    model_handles = [Line2D([0], [0], **PL_STYLE, label=r"\makebox[2.3em][l]{$F^{\mathrm{PL}}$}point-like"),
                     Line2D([0], [0], color="0.2", **WS_STYLE, label=r"\makebox[2.3em][l]{$F$}at $\mathbf{b}$"),
                     Line2D([0], [0], color="0.2", **EFF_STYLE, label=r"\makebox[2.3em][l]{$F^{\mathrm{eff}}$}at $\mathbf{r}$")]
    for a in (ax, ax_ratio):
        a.tick_params(labelsize=26)
        a.tick_params(axis="x", pad=9)
        a.tick_params(axis="y", pad=6)
    channel_legend = ax.legend(handles=channel_handles, loc="lower left", bbox_to_anchor=(0.04, 0.0), fontsize=22, labelspacing=0.45, frameon=False)
    ax.add_artist(channel_legend)
    ax.legend(handles=model_handles, loc="lower left", bbox_to_anchor=(0.04, 0.325), fontsize=22, labelspacing=0.45, frameon=False)

    # lower panel: the three F^eff/F lines lie on top of each other, so list them in one row
    ratio_handles = [Line2D([0], [0], color=CHANNEL_COLORS[ch], **EFF_STYLE, label=ch) for ch in CHANNELS]
    ratio_handles.append(Line2D([0], [0], **PL_STYLE, label=r"$F^{\mathrm{PL}}$"))
    ax_ratio.legend(handles=ratio_handles, loc="upper left", bbox_to_anchor=(0.0, 0.88), ncol=4, fontsize=22, frameon=False,
                    handlelength=1.0, handletextpad=0.3, columnspacing=0.6, borderaxespad=0.3)
    fig.align_ylabels([ax, ax_ratio])   # both y titles at the same distance from the axis
    outname = "../plots/flux_eff_vs_ws.pdf"
    plt.savefig(outname, bbox_inches="tight")
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
