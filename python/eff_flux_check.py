import csv, os, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

plt.rcParams.update({
    "text.usetex": True, "font.family": "serif", "font.size": 11,
    "axes.labelsize": 20, "xtick.labelsize": 18, "ytick.labelsize": 18,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 8, "ytick.major.size": 8,
    "xtick.minor.size": 4, "ytick.minor.size": 4,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "xtick.color": "0.4", "ytick.color": "0.4",
    "xtick.major.pad": 8, "ytick.major.pad": 4,
})

def load_table(path):
    rows = [l.split() for l in open(path) if l.strip() and not l.startswith("#")]
    a = np.array(rows, dtype=float)
    return {"y4": a[:,2], "AnAn": a[:,3], "An0n": a[:,4]}

def interp_flux(tab, key, y):
    y4j, logf = tab["y4"], tab[key]
    n = len(y4j)-1
    w = (-1.0)**np.arange(n+1); w[0]*=0.5; w[-1]*=0.5
    out=[]
    for yy in np.atleast_1d(y):
        d = yy**0.25 - y4j
        k = np.argmin(abs(d))
        if abs(d[k])<1e-14: out.append(np.exp(logf[k])); continue
        r = w/d
        out.append(np.exp(np.sum(r*logf)/np.sum(r)))
    return np.array(out)

def load_scan(fn):
    d={}
    for r in csv.DictReader(open(fn)):
        d.setdefault(r["channel"],[]).append((float(r["y"]),float(r["dN_dy"])))
    return {k:(np.array([p[0] for p in v]), np.array([p[1] for p in v])) for k,v in d.items()}

WS = load_table("../inputs/photon_flux/log-flux-tbl-WS.dta")
simple = load_scan("../out/flux_scan_PL.csv")
eff = {}
eff["AnAn"] = load_scan("../out/flux_scan_eff_AnAn.csv")["AnAn"]
eff["An0n"] = load_scan("../out/flux_scan_eff_An0n.csv")["An0n"]

fig, ax = plt.subplots(figsize=(7.5,6))
yy = np.logspace(-4, -0.3, 300)
colors = {"AnAn": "#2a78d6", "An0n": "#eb6834"}
for ch in ["AnAn", "An0n"]:
    y_s, f_s = simple[ch]
    simple_ratio = np.exp(np.interp(np.log(yy), np.log(y_s), np.log(f_s))) / interp_flux(WS, ch, yy)
    y_e, f_e = eff[ch]
    eff_ratio = np.exp(np.interp(np.log(yy), np.log(y_e), np.log(f_e))) / interp_flux(WS, ch, yy)
    ax.plot(yy, simple_ratio, color=colors[ch], ls=":", lw=2, label=f"{ch}, simple (Eq. 1)")
    ax.plot(yy, eff_ratio, color=colors[ch], ls="-", lw=2, label=f"{ch}, effective (Eq. 4)")

ax.axhline(1.0, color="0.6", lw=1, ls="--")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1e-4, 0.5); ax.set_ylim(1e-2, 3)
ax.set_xlabel(r"$y=2\omega/\sqrt{s_{NN}}$", labelpad=15)
ax.set_ylabel("this code / Starlight table", labelpad=15)
for a in [ax]:
    a.xaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    a.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2,10)*0.1, numticks=20))
    a.yaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    a.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2,10)*0.1, numticks=20))
ax.legend(fontsize=12, frameon=False, loc="lower left")
ax.text(0.95,0.95,"Pb-Pb 5.36 TeV", transform=ax.transAxes, ha="right", va="top", fontsize=14)
plt.tight_layout()
plt.savefig("../plots/eff_flux_check.pdf")
print("Saved ../plots/eff_flux_check.pdf")
