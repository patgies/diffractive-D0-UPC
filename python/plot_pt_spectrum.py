import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from alphas_running import alphas_run

# Plots the D0-level dsigma/dpD0 spectrum vs pD0 at each fixed rapidity y,
# reading the flat (pD0, y, dsigma) files written by ../run_local.sh
# (src/main.cpp / d0_point, BCFY fragmentation + An0n channel):
#   ../files/d0_point_exclusive.dat
#   ../files/d0_point_diffractive.dat

alphae = 1/137
mc     = 1.5          # charm mass GeV
e_c    = 2/3
Nc     = 3
sigma0 = 16.36      # mb, dipole normalization
Y_FLOOR = 1e-12


def load_points(fname):
    """Returns (pD0, y, dsigma) arrays, one row per computed point."""
    data = np.loadtxt(fname, comments='#')
    return data[:, 0], data[:, 1], data[:, 2]


excl_file = "../files/d0_point_exclusive.dat"
diff_file = "../files/d0_point_diffractive.dat"

if not os.path.exists(excl_file) or not os.path.exists(diff_file):
    sys.exit(f"Missing file(s): {excl_file}, {diff_file} "
             "(run ../run_local.sh first)")

pD0_excl, y_excl, excl = load_points(excl_file)
pD0_diff, y_diff, diff = load_points(diff_file)
assert np.allclose(pD0_excl, pD0_diff) and np.allclose(y_excl, y_diff), \
    "exclusive/diffractive points don't match"

y_values = np.unique(y_excl)
colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(y_values)))

fig, ax = plt.subplots(figsize=(7.5, 6.5))

# Same missing-prefactor convention as the parton-level integrands: the
# fragmented cross sections returned by d0_point still need
# alpha_s * alpha_em * e_c^2 * (Nc^2-1) * sigma0 / (8*pi^4) (diffractive) and
# alpha_em * Nc * e_c^2 * sigma0 / (2*pi^2) (exclusive) applied downstream.
p_e = alphae * Nc * e_c**2 * sigma0 / (2 * np.pi**2)

for y, c in zip(y_values, colors):
    mask = np.isclose(y_excl, y)
    order = np.argsort(pD0_excl[mask])
    pt = pD0_excl[mask][order]

    alphas = alphas_run(np.sqrt(pt**2 + mc**2))
    p_d = alphas * alphae * e_c**2 * (Nc**2 - 1) * sigma0 / (8 * np.pi**4)

    dsig_diff = 2*np.pi * pt * diff[mask][order] * p_d
    dsig_excl = 2*np.pi * pt * excl[mask][order] * p_e

    ax.plot(pt, dsig_diff, color=c, linestyle='-',  label=f"diff., $y={y:g}$")
    ax.plot(pt, dsig_excl, color=c, linestyle='--', label=f"excl., $y={y:g}$")

ax.set_yscale('log')
ax.set_xlabel(r'$p_{D^0}$ [GeV]')
ax.set_ylabel(r'$d\sigma/dp_{D^0}$ (mb/GeV, arb. norm.)')
ax.set_title(r'D0-level (fragmented) $p_T$ spectrum')
ax.tick_params(axis='both', which='major', labelsize=12)

# Reorder so the diffractive entries fill the left column and the
# exclusive entries fill the right one (see plot_flux_spectrum.py).
handles, labels = ax.get_legend_handles_labels()
diff_items = [(h, l) for h, l in zip(handles, labels) if l.startswith('diff')]
excl_items = [(h, l) for h, l in zip(handles, labels) if l.startswith('excl')]
handles, labels = zip(*(diff_items + excl_items))
ax.legend(handles, labels, fontsize=7.5, ncol=2)
ax.grid(True, which='both', alpha=0.3)
ax.set_ylim(bottom=Y_FLOOR)

plt.tight_layout()
outname = "../plots/pt_spectrum.png"
plt.savefig(outname, dpi=300)
print(f"Saved: {outname}")
