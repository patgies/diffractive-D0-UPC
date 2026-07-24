import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from alphas_running import alphas_run

# Plots dsigma/dpt vs pt at fixed rapidity y, reading the 2D (pt,y) grid
# files written by ../build/bin/flux_grid (main_grid.cpp):
#   ../files/flux_grid_exclusive.dat
#   ../files/flux_grid_diffractive.dat

alphae = 1/137
mc     = 1.5          # charm mass GeV
e_c    = 2/3
Nc     = 3
sigma0 = 16.36      # mb, dipole normalization
Y_FLOOR = 1e-12       


def load_grid(fname):
    """Returns (pt, y_grid, values[len(pt), len(y_grid)])."""
    y_grid = None
    with open(fname) as fh:
        for line in fh:
            if line.startswith('# pt \\ y'):
                y_grid = np.array([float(v) for v in line.split()[4:]])
                break
    if y_grid is None:
        raise ValueError(f"Could not find 'pt \\ y' header line in {fname}")
    data = np.loadtxt(fname, comments='#')
    pt = data[:, 0]
    values = data[:, 1:]
    return pt, y_grid, values


excl_file = "../files/flux_grid_exclusive.dat"
diff_file = "../files/flux_grid_diffractive.dat"

if not os.path.exists(excl_file) or not os.path.exists(diff_file):
    sys.exit(f"Missing grid file(s): {excl_file}, {diff_file} "
             "(run ./run_local.sh first)")

pt, y_excl, excl = load_grid(excl_file)
_, y_diff, diff = load_grid(diff_file)
assert np.allclose(y_excl, y_diff), "exclusive/diffractive y grids don't match"
y_values = y_excl

colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(y_values)))

fig, ax = plt.subplots(figsize=(7.5, 6.5))

alphas = alphas_run(np.sqrt(pt**2 + mc**2))
p_d = alphas * alphae * e_c**2 * (Nc**2 - 1) * sigma0 / (8 * np.pi**4)
p_e = alphae * Nc * e_c**2 * sigma0 / (2 * np.pi**2)

for j, (y, c) in enumerate(zip(y_values, colors)):
    dsig_diff = 2*np.pi * pt * diff[:, j] * p_d
    dsig_excl = 2*np.pi * pt * excl[:, j] * p_e

    ax.plot(pt, dsig_diff, color=c, linestyle='-',  label=f"diff., $y={y:g}$")
    ax.plot(pt, dsig_excl, color=c, linestyle='--', label=f"excl., $y={y:g}$")

ax.set_yscale('log')
ax.set_xlabel(r'$p_T$ [GeV]')
ax.set_ylabel(r'$d\sigma/dp_T$ (mb/GeV, arb. norm.)')
ax.set_title(r'Photon-flux-convolved $p_T$ spectrum (no fixed $q^+$)')
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
outname = "../plots/flux_grid_spectrum.png"
plt.savefig(outname, dpi=300)
print(f"Saved: {outname}")
