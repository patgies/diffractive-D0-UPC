import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Same "fancy" styling as cross_section.py / xpom_spectrum.py.
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
})

mc = 1.5   # charm mass in GeV, same as cross_section.py's convention

pt_values = np.linspace(0.0, 8.0, 400)
mt_values = np.sqrt(mc**2 + pt_values**2)

plt.figure(figsize=(7.5, 6.0))
plt.plot(pt_values, mt_values, color="#08519c", linewidth=2, label=r"$m_T=\sqrt{m_c^2+p_{D^0\perp}^2}$")
plt.plot(pt_values, pt_values, color="gray", linestyle="--", linewidth=1.5, label=r"$m_T=p_{D^0\perp}$ (asymptote)")
plt.axhline(mc, color="gray", linestyle=":", linewidth=1)
plt.text(6.5, mc + 0.15, r"$m_c$", fontsize=16, color="gray")

plt.xlabel(r"$p_{D^0\perp}$ [GeV]", labelpad=15)
plt.ylabel(r"$m_T$ [GeV]", labelpad=15)
plt.title(r"Transverse mass $m_T=\sqrt{m_c^2+p_{D^0\perp}^2}$", pad=15)
plt.legend(loc="upper left", fontsize=15)

plt.figtext(0.5, -0.02,
            r"At low $p_{D^0\perp}$, $m_T\to m_c=1.5$~GeV (its minimum); at high $p_{D^0\perp}$, $m_T\to p_{D^0\perp}$.",
            fontsize=9, color="dimgray", ha="center")

plt.tight_layout()
outname = "../plots/mt_vs_pt.pdf"
plt.savefig(outname, bbox_inches="tight")
print(f"Saved: {outname}")
