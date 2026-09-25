import csv
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Compares the b-integrated photon flux of this code (src/scan_flux.cpp ->
# out/flux_scan.csv, dN/dy with y = 2*omega/sqrt(s_NN)) with the tabulated
# effective fluxes of arXiv:2404.09731 (inputs/photon_flux/log-flux-tbl-*.dta,
# P. Paakkinen), evaluated with the paper's barycentric Lagrange
# interpolation on Chebyshev nodes in y^(1/4) (Eqs. 34-37).
# Run from python/:  ../build/bin/scan_flux > ../out/flux_scan.csv first.

TABLE_DIR = "../inputs/photon_flux"
SCAN = "../out/flux_scan.csv"


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


def load_scan():
    d = {}
    for r in csv.DictReader(open(SCAN)):
        d.setdefault((r["channel"], float(r["bmax"])), []).append((float(r["y"]), float(r["dN_dy"])))
    return {k: (np.array([p[0] for p in v]), np.array([p[1] for p in v])) for k, v in d.items()}


PL = load_table(os.path.join(TABLE_DIR, "log-flux-tbl-PL.dta"))
WS = load_table(os.path.join(TABLE_DIR, "log-flux-tbl-WS.dta"))
scan = load_scan()

# (our channel, table, table key, label, color)
PAIRS = [
    ("PL(AnAn)", PL, "AnAn", "PL, AnAn", "tab:blue"),
    ("AnAn",     WS, "AnAn", r"$\Gamma_{AA}$ survival vs WS, AnAn", "tab:green"),
    ("An0n",     WS, "An0n", r"$\Gamma_{AA}$ + EMD vs WS, An0n", "tab:red"),
]

if __name__ == "__main__":
    print(f"{'y':>9} " + " ".join(f"{c[3][:14]:>16}" for c in PAIRS) + "   (ours[bmax=1e7]/table)")
    for yy in [1e-4, 1e-3, 1e-2, 1e-1, 0.5]:
        line = f"{yy:9.1e} "
        for ch, tab, key, _, _ in PAIRS:
            y, f = scan[(ch, 1.0e7)]
            ours = np.interp(np.log(yy), np.log(y), np.log(f))
            line += f"{np.exp(ours) / interp_flux(tab, key, yy)[0]:16.3f} "
        print(line)

    fig, (ax, axr) = plt.subplots(2, 1, figsize=(7.5, 8), sharex=True, gridspec_kw={"height_ratios": [3, 2]})
    yy = np.logspace(-4, 0, 400)
    for ch, tab, key, label, color in PAIRS:
        ax.plot(yy, yy * interp_flux(tab, key, yy), color=color, lw=2, label=f"table: {label}")
        for bmax, ls in [(1.0e7, "--"), (650.0, ":")]:
            y, f = scan[(ch, bmax)]
            ax.plot(y, y * f, color="k" if False else color, ls=ls, lw=1.3,
                    label=f"this code, b<{'inf' if bmax > 1e6 else '650 GeV$^{-1}$'}" if ch == "PL(AnAn)" else None)
            axr.plot(y, f / interp_flux(tab, key, y), color=color, ls=ls, lw=1.5)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylim(1e-9, 1e3); ax.set_xlim(1e-4, 0.5)
    ax.set_ylabel(r"$y\,f(y) = dN/d\ln y$")
    ax.set_title(r"Pb+Pb $\sqrt{s_{NN}}=5.36$ TeV: $b$-integrated photon flux")
    ax.legend(fontsize=8)
    axr.axhline(1.0, color="gray", lw=1)
    axr.set_xscale("log"); axr.set_ylabel("this code / table"); axr.set_xlabel(r"$y=2\omega/\sqrt{s_{NN}}$")
    axr.set_yscale("log"); axr.set_ylim(1e-3, 3)
    plt.tight_layout()
    out = "../plots/flux_comparison.pdf"
    plt.savefig(out)
    print("Saved:", out)
