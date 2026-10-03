# Diffractive D0 photoproduction

This project calculates the diffractive D0 photoproduction cross section `dσ / (dy d²p_D0)` in ultraperipheral collisions (UPCs) in the CGC framework. The code supports different UPC channels, photon fluxes and fragmentation functions. It can be used for proton and nuclear targets.

If you use this code, please cite:

P. Caucal, P. Gimeno-Estivill, E. Iancu, T. Lappi, and F. Salazar, *TMD factorization in diffractive heavy quark production in photon-nucleus collisions*, JHEP 09 (2026) 169 [[arXiv:2606.04169](https://arxiv.org/abs/2606.04169)] [[doi:10.1007/JHEP09(2026)169](https://doi.org/10.1007/JHEP09(2026)169)].

---

## Build

```bash
mkdir build
cd build
cmake ..
make
```

Requirements:

- CMake
- GSL (GNU Scientific Library)
- OpenMP

---

## Basic run

The C++ sources are in [src](src), the run scripts in [run_scripts](run_scripts) and the plotting scripts in [plotting_scripts](plotting_scripts). The main executable is

```bash
./build/bin/D0 <dipole_file> <pD0> <y> [<FRAG_TYPE>] [<CHANNEL>]
```

Example:

```bash
./build/bin/D0 data/proton/mve.dat 3 0 HymnD An0n
```

This prints `pD0`, the exclusive and the diffractive cross section, without prefactors.

The other executables are `D0_xpom` (differential in `x_pom`), `charm_fixed_qp` (charm quark at fixed `q+`, no photon flux) and `scan_flux` (photon flux alone).

---

## Run scripts

The results go to `output/`.

- `run_nucleus.sh`: nuclear target, every `b_d`, `pD0` and `y`
- `run_scale_variation.sh`: the same with the fragmentation scale `Q = 0.5 mT` and `2 mT`
- `run_proton.sh`: `TARGET=pA` for p+Pb at 8.16 TeV
- `run_xpom.sh`: cross section differential in `x_pom`
- `run_members.sh`: HymnD replicas (`MEMBER_SET=HymnD`) or BK posterior samples (`MEMBER_SET=bk`)
- `run_charm_fixed_qp.sh`: charm quark at fixed `q+`
- `run_flux_scan.sh`: photon flux alone

The settings are environment variables, defaults are in [run_scripts/config.sh](run_scripts/config.sh).

```bash
CHANNEL=0n0n FRAG_TYPE=BCFY ./run_scripts/run_nucleus.sh
Y_VALS="0.0 1.0" PT_VALS="1.0 2.0 4.0" CALLS=1e4 ./run_scripts/run_nucleus.sh
```

---

## Inputs and targets

The code expects dipole input files of the form

- proton: `data/proton/mve.dat`
- nucleus: `data/Pb/mve/glauber_mve_<b_d>` or `data/Au/...`, one Glauber-sampled file per impact parameter `b_d`

Dipole parametrization MVe from [https://github.com/hejajama/rcbkdipole](https://github.com/hejajama/rcbkdipole).

The other inputs are in `input/`:

- `BCFY_EKO/` and `KK_EKO/`: grids of the BCFY and Kniehl-Kramer fragmentation functions, DGLAP-evolved with [eko](https://github.com/NNPDF/eko) in [EKO-FF](https://github.com/patgies/EKO-FF).
- `HymnD/`: the HymnD fragmentation function set, members 0–100.
- `WS_photon_flux/`: tables of the hadronic survival factor `Gamma_AA(b)`, for `sigma_NN = 92` mb (default) and `90.85` mb. They are made with [src/make_gamma_aa.py](src/make_gamma_aa.py).
- `Starlight_photon_flux/`: the Starlight photon flux tables, used with `FLUX_MODEL=STARLIGHT` and in `flux_comparison.py`.
- `BK/`: the posterior samples of the BK initial condition. `run_scripts/setup_bk_posterior_links.sh` links them as `data/Pb/bk_posterior/member_NNNN/`.
- `inclusive/`: the inclusive D0 cross section from [inclusive-D0-UPC](https://github.com/patgies/inclusive-D0-UPC), only needed for `ratio.py`. It is not included.

---

## UPC channel and photon flux

The channel is selected through `CHANNEL`: `An0n`, `Xn0n`, `0n0n`, `AnAn`, `PL(AnAn)`.

The photon flux is selected through `FLUX_MODEL`:

- `EFF` (default): effective flux of K. J. Eskola, V. Guzey, I. Helenius, P. Paakkinen, and H. Paukkunen, "Spatial resolution of dijet photoproduction in near-encounter ultraperipheral nuclear collisions," Phys. Rev. C 110, 054906 (2024) [arXiv:2404.09731], Eq. (4).
- `STARLIGHT`: flux tables of the same paper (only `AnAn` and `An0n`).
- `PL`: flux of a point-like nucleus.
- `WS`: flux for a Woods-Saxon charge distribution.

The photon flux and how to change it is explained in [notes_photon_flux.pdf](notes_photon_flux.pdf).

---

## Fragmentation functions

The fragmentation function is selected through `FRAG_TYPE`:

- `BCFY`: E. Braaten, K.-m. Cheung, S. Fleming, and T.-C. Yuan, "Perturbative QCD fragmentation functions as a model for heavy quark fragmentation," Phys. Rev. D 51, 4819 (1995) [arXiv:hep-ph/9409316].
- `KniehlKramer`: B. A. Kniehl and G. Kramer, "Charmed-hadron fragmentation functions from CERN LEP1 revisited," Phys. Rev. D 74, 037502 (2006) [arXiv:hep-ph/0607306].
- `HymnD` (default): Epele, Hekhorn, Helenius, Paukkunen, and Zurita, "Towards new D meson fragmentation functions" [arXiv:2609.10327].

`BCFY` and `KniehlKramer` are evolved with DGLAP from their input at the starting scale `mc` up to the fragmentation scale. The evolution is done with [eko](https://github.com/NNPDF/eko) in the repository [EKO-FF](https://github.com/patgies/EKO-FF) and stored as grids in `input/BCFY_EKO/` and `input/KK_EKO/`.

The fragmentation scale is `Q = SCALE_FACTOR * mT`.

The fragmentation functions and how to change them is explained in [notes_FF.pdf](notes_FF.pdf).

The HymnD members 0–100 are in `input/HymnD/`.

### Uncertainty bands

The plots include the scale uncertainty for the three fragmentation functions (`HymnD`, `BCFY` and `KniehlKramer`). [D0_bins.py](plotting_scripts/D0_bins.py) shows it for all three; [D0.py](plotting_scripts/D0.py) and [D0_sum.py](plotting_scripts/D0_sum.py) show it only for the HymnD curve:

- a band from varying the fragmentation scale `Q` and the renormalization scale of `alpha_s` by a factor of `0.5` and `2` around the central scale `mT`. The band is the envelope of the 7 combinations where the two scales differ by at most a factor of 2.

The fragmentation scale variation is run with `run_scale_variation.sh`, once per fragmentation function. The renormalization scale variation is done in the plotting scripts and only changes the diffractive cross section.

Two other uncertainties can be computed with `run_members.sh` but are not included in the plots:

- the HymnD replicas (`MEMBER_SET=HymnD`): one run per member of the fragmentation function set.
- the BK initial condition (`MEMBER_SET=bk`): one run per posterior sample of the BK fit, with the fragmentation function fixed. The samples are in `input/BK/` and `data/Pb/bk_posterior/`.

---

## Differential cross section and normalization

The output of the nuclear runs is per impact parameter `b_d`. The plotting scripts integrate over `b_d` (Simpson's rule, multiplied by `2πb_d`) and multiply by the prefactors.

For a proton target, the impact-parameter integral is replaced by the proton normalization `16.36` mb of the MVe dipole parametrization.

Units are GeV throughout.
