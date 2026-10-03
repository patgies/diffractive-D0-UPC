# diffractive-D0-UPC

Exclusive and diffractive D⁰ photoproduction in ultraperipheral Pb+Pb collisions at
√s<sub>NN</sub> = 5.36 TeV, in the color-dipole picture. The C++ programs compute
dσ/(dy d²p<sub>D⁰⊥</sub>) for a given D⁰ transverse momentum and rapidity; the
shell workflows scan these over a (p<sub>T</sub>, y) grid and the impact parameter
b<sub>d</sub> of the dipole in the target (one Glauber sample per value), and the Python
scripts integrate over b<sub>d</sub> and make the plots.

Included:

- **Photon flux:** effective flux with hadronic survival factor Γ<sub>AA</sub>(b) (σ<sub>NN</sub> = 92 mb)
  and the neutron classes **An0n**, **Xn0n**, **0n0n** and AnAn (electromagnetic dissociation
  probability P<sub>noEM</sub>(b) = exp(−S/b²), squared for 0n0n); a point-like flux and the
  Starlight tables for comparison.
- **Fragmentation:** BCFY, Kniehl–Kramer and HymnD fragmentation functions, with the
  fragmentation-scale variation Q = {0.5, 1, 2} m<sub>T</sub>.
- **Uncertainty bands:** scale variation, HymnD replicas and BK initial-condition posterior samples.
- **Other observables:** the R<sub>pA</sub> proton baseline, p+Pb at 8.16 TeV, x<sub>ℙ</sub> spectra and
  a fixed-q⁺ charm-quark cross section without photon flux.

## Notation

The code uses the same names as the notes:

| Name in the code | | Meaning |
|---|---|---|
| `z`, `zbar` | z, 1 − z | momentum fraction of the photon carried by the charm quark, z = p⁺/q⁺ |
| `z_gamma` | z<sub>γ</sub> | energy fraction of the beam nucleon carried by the photon, z<sub>γ</sub> = 2ω/√s<sub>NN</sub> (the Starlight tables call it y) |
| `z_h` | z<sub>h</sub> | fragmentation fraction, p<sub>D⁰</sub> = z<sub>h</sub> p<sub>c</sub>, integrated from `z_h_min` to `z_h_max` |
| `b_d` | b<sub>d</sub> | impact parameter of the dipole in the target. Not a C++ variable: it is fixed by the Glauber file `glauber_mve_<b_d>` and integrated by the Python scripts |
| `b` | b | distance between the centres of the two nuclei; Γ<sub>AA</sub>(b) and the neutron-class factors depend on it |
| `r` | r | distance from the photon to the centre of the nucleus that emits it (argument of the flux densities) |
| `s` | s | position in the target nucleus, b = \|r − s\| (only in the effective flux) |

## Requirements

- C++14 compiler, CMake ≥ 3.0, [GSL](https://www.gnu.org/software/gsl/), OpenMP
- bash, awk, `seq` (for the workflows)
- Python 3 with numpy, scipy and matplotlib, and a LaTeX installation (the plots use `text.usetex`)

## Build

```bash
cmake -S . -B build
cmake --build build -j
```

This builds `build/bin/D0` (the main program), `D0_xpom` (x<sub>ℙ</sub>-differential),
`charm_fixed_qp`, and the flux check `scan_flux`. The workflow
scripts do not build anything: build first, and again after changing the code.

## Running the calculations

All workflows are in `run_scripts/` and run on an ordinary machine (no cluster
needed). They can be started from any folder; results go to `output/`.

Each calculation has its own script:

| Step | Script | Output |
|---|---|---|
| `central` | `run_nucleus.sh` | `output/<CHANNEL>/central_values/<FRAG>/` |
| `scale` | `run_scale_variation.sh` | `output/<CHANNEL>/scale_variation/<FRAG>/factor_{0.5,2.0}/` |
| `proton` | `run_proton.sh` | `output/<CHANNEL>/proton_baseline/<FRAG>/` (R<sub>pA</sub> baseline) |
| `xpom` | `run_xpom.sh` | `output/<CHANNEL>/xpom/<FRAG>/` |
| `hymnd_band` | `run_members.sh` (`MEMBER_SET=HymnD`) | `output/<CHANNEL>/HymnD_band/member_NNNN/` |
| `bk_band` | `run_members.sh` (`MEMBER_SET=bk`) | `output/<CHANNEL>/bk_band/member_NNNN/` |
| `pPb` | `run_proton.sh` (`TARGET=pA`) | `output/pPb/<FRAG>/` (p+Pb at 8.16 TeV, one run per FF in `FRAGS`) |
| `charm` | `run_charm_fixed_qp.sh` | `output/charm/` |
| `flux` | `run_flux_scan.sh` | `output/flux_scan/` (photon flux alone, for `flux_comparison.py`; σ<sub>NN</sub> = 90.85 mb of the Starlight tables, change it with `GAMMA_AA_FILE`) |

Every result is in a subfolder named after its fragmentation function (`<FRAG>`), except
`charm/` and `flux_scan/`, which use none. Run a script with its settings on the command line, e.g.
`CHANNEL=0n0n FRAG_TYPE=BCFY ./run_scripts/run_nucleus.sh`.
`run_proton.sh` with `TARGET=pA` (the `pPb` step) computes p+Pb collisions at 8.16 TeV.

### Settings

Every setting is an environment variable with its default in
[`run_scripts/config.sh`](run_scripts/config.sh); set it on the command line to
change it. The main ones:

| Variable | Default | |
|---|---|---|
| `CHANNEL` / `CHANNELS` | `An0n` | neutron class: `An0n`, `Xn0n`, `0n0n`, `AnAn`, `PL(AnAn)` |
| `FRAG_TYPE` / `FRAGS` | `HymnD` / all three | `BCFY`, `KniehlKramer`, `HymnD` |
| `Y_VALS` | −2.0 … 2.0, step 0.5 | rapidities |
| `PT_VALS` | 0.2 … 2.0 (step 0.1), 2.5 … 12.0 (step 0.5) | p<sub>D⁰⊥</sub> in GeV |
| `FLUX_MODEL` | `EFF` | `EFF`, `STARLIGHT`, `PL`, `WS` |
| `CALLS`, `CALLS_EXCL_FACTORIZED` | `1e5`, `2e3` | VEGAS calls per point |
| `CORES` | half the CPUs | parallel processes |
| `OUTPUT_ROOT` | `output` | folder all results go under |

### Using a different photon flux

The flux is chosen with `FLUX_MODEL` (`EFF`, `STARLIGHT`, `PL`, `WS`), the neutron class
with `CHANNEL`, and σ<sub>NN</sub> through the Γ<sub>AA</sub> table in `GAMMA_AA_FILE` (default
92 mb; `input/WS_photon_flux/Gamma_AA_sigma90.85.dat` has the 90.85 mb of the Starlight tables).

To add your own flux, follow one rule: **every flux function in `src/photon_flux.cpp`
returns dN/dω**. Convert your table there (e.g. multiply dN/dz<sub>γ</sub> by 2/√s<sub>NN</sub>), and
never add flux-specific factors of z<sub>γ</sub> or 1/q⁺ in the integrands. Then check it with
`./build/bin/scan_flux <FLUX_MODEL> <CHANNEL>`, which must reproduce your table, and by
comparing one D⁰ point with two fluxes. Step-by-step instructions, the formulas and the
checks are in [`notes_photon_flux.pdf`](notes_photon_flux.pdf).

### Run time

A full central run (one channel, one fragmentation function, the full grid) is about
6000 program calls over the 17 b<sub>d</sub> samples: roughly 30 min on 64 cores,
or several hours on a laptop. To try things out, use a smaller grid and fewer calls:

```bash
Y_VALS="0.0 1.0" PT_VALS="1.0 2.0 4.0" CALLS=1e4 ./run_scripts/run_nucleus.sh
```

`run_proton.sh` and `run_charm_fixed_qp.sh` are cheap (minutes). The member bands are the
most expensive: each member is a full central run.

On a cluster, each step can be run as a separate job: the scripts take all their
settings from the environment, so a job only needs to export them and call the script
(with `CORES` set to the job's CPU count).

## Plots

Run the plotting scripts from `plotting_scripts/`; they read `output/<CHANNEL>/` and write to
`plots/`. Use `CHANNEL=...` for another neutron class.

| Script | Plot | Needs |
|---|---|---|
| `D0.py` | `D0_<CHANNEL>_PbPb.pdf`: exclusive and diffractive spectra with scale band | `central`, `scale` |
| `D0_sum.py` | `D0_sum_<CHANNEL>_PbPb.pdf`: exclusive + diffractive | `central`, `scale` |
| `D0_bins.py` | `D0_bins_<CHANNEL>_PbPb.pdf`: averages over (p<sub>T</sub>, y) bins | `central`, `scale` (all FFs) |
| `RpA.py` | `RpA_<CHANNEL>.pdf`, `RpA_<CHANNEL>_exclusive.pdf` | `central`, `proton` |
| `xpom.py` | `xpom_<CHANNEL>_PbPb.pdf`: diffractive x<sub>ℙ</sub> spectrum | `xpom` |
| `ratio.py` | `ratio_<CHANNEL>_PbPb.pdf`: (diffractive + exclusive) / inclusive | `central` (BCFY) and `input/inclusive/` |
| `charm_fixed_qp.py` | `charm_fixed_qp.pdf` | `charm` |
| `flux_comparison.py` | `flux_comparison.pdf`: photon fluxes vs. the Starlight tables | `output/flux_scan/` |

```bash
cd plotting_scripts
python3 D0.py
CHANNEL=0n0n python3 D0_sum.py
```

`src/combine_nucleus.py` (run it from `src/`) writes the b<sub>d</sub>-integrated spectra as two-column
tables in `output/<CHANNEL>/combined/`.

## Repository layout

| Folder | Contents |
|---|---|
| `src/` | C++ source; `make_gamma_aa.py` (makes the Γ<sub>AA</sub>(b) tables in `input/WS_photon_flux/`), `combine_nucleus.py` (spectra as tables) and `alphas_running.py` (running coupling, used by the plotting scripts) |
| `run_scripts/` | the shell workflows (see above) |
| `plotting_scripts/` | b<sub>d</sub> integration and plots |
| `data/` | dipole amplitudes: Glauber samples for Pb and Au (`<NUCLEUS>/mve/glauber_mve_<b_d>`), the proton, and the BK posterior samples |
| `input/` | fragmentation-function grids, photon-flux tables, BK posterior inputs |
| `output/` | results of the workflows (the current results are included) |
| `plots/` | plots |
| `notes/` | notes on the formulae, the photon flux and the fragmentation functions (`.tex` + `.pdf`) |
| `notes_photon_flux.pdf` | the photon flux, and how to use your own (source in `notes/`) |
| `notes_FF.pdf` | the fragmentation functions, and how to choose one or add your own (source in `notes/`) |

## Not included

- **HymnD replicas.** Only the central member `input/HymnD/prompt-D0-1-109_0000.dat` is
  included. For the HymnD replica band, place the members 1–100 as
  `input/HymnD/prompt-D0-1-109_NNNN.dat`.
- **Inclusive D⁰ cross section.** `ratio.py` needs the inclusive cross section from the
  separate inclusive-D0-UPC project, in `input/inclusive/`.

The BK posterior samples are included (`input/BK/`); the links in
`data/Pb/bk_posterior/` can be recreated with `run_scripts/setup_bk_posterior_links.sh`.
