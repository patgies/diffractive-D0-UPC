import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from alphas_running import alphas_run, mZ

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

mu_values = np.logspace(np.log10(mc), np.log10(mZ), 400)
alphas_values = alphas_run(mu_values)

plt.figure(figsize=(7.5, 6.0))
plt.plot(mu_values, alphas_values, color="#08519c", linewidth=2)

# Mark the charm-mass threshold and the m_Z normalization point used to fix
# the one-loop running (alpha_s(mZ)=0.118, PDG world average).
plt.axvline(mc, color="gray", linestyle=":", linewidth=1)
plt.text(mc * 1.15, 0.30, r"$m_c$", fontsize=16, color="gray")
plt.axvline(mZ, color="gray", linestyle=":", linewidth=1)
plt.text(mZ * 0.7, 0.30, r"$m_Z$", fontsize=16, color="gray")

# The p_D0 range used elsewhere in this project maps to mu=sqrt(mc^2+pD0^2);
# shade that range for reference (pD0 up to 8 GeV, as in cross_section.py).
pt_max = 8.0
mu_max_used = np.sqrt(mc**2 + pt_max**2)
plt.axvspan(mc, mu_max_used, color="#9ecae1", alpha=0.2, linewidth=0)

plt.xscale("log")
plt.xlabel(r"$\mu$ [GeV]", labelpad=15)
plt.ylabel(r"$\alpha_s(\mu)$", labelpad=15)
plt.title(r"One-loop running coupling ($\alpha_s(m_Z)=0.118$, $N_f=4$)", pad=15)

plt.figtext(0.5, -0.02,
            r"Shaded: $\mu=\sqrt{m_c^2+p_{D^0\perp}^2}$ range for $p_{D^0\perp}\in[0,8]$~GeV, "
            r"the scale used in this project's cross sections.",
            fontsize=9, color="dimgray", ha="center")

plt.tight_layout()
outname = "../plots/alphas_running.pdf"
plt.savefig(outname, bbox_inches="tight")
print(f"Saved: {outname}")
