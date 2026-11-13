# Diffractive D0 photoproduction

Computes the exclusive and diffractive D0 photoproduction cross section `dσ / (dy d²p_D0)` in ultraperipheral collisions (UPCs), in the color glass condensate framework. It works for proton and nuclear targets, with several UPC channels, photon fluxes and fragmentation functions.

If you use this code, please cite:

P. Caucal, P. Gimeno-Estivill, E. Iancu, T. Lappi, and F. Salazar, *TMD factorization in diffractive heavy quark production in photon-nucleus collisions*, JHEP 09 (2026) 169 [[arXiv:2606.04169](https://arxiv.org/abs/2606.04169)] [[doi:10.1007/JHEP09(2026)169](https://doi.org/10.1007/JHEP09(2026)169)].

---

## Build

```bash
mkdir build && cd build
cmake ..
make
```

Requirements: CMake, GSL (GNU Scientific Library) and OpenMP.

## Repository layout

- [src](src): C++ sources
- [run_scripts](run_scripts): run scripts, with the default settings in [config.sh](run_scripts/config.sh)
- [plotting_scripts](plotting_scripts): plotting scripts
- `data/`: dipole files, `input/`: other inputs, `output/`: results, `plots/`: figures

## Basic run

```bash
./build/bin/D0 <dipole_file> <pD0> <y> [<FRAG_TYPE>] [<CHANNEL>]
./build/bin/D0 data/proton/mve.dat 3 0 HymnD An0n
```

This prints `pD0`, the exclusive and the diffractive cross section, without prefactors.

Other executables:

| Executable | Computes |
|---|---|
| `D0_xpom` | the D0 cross section, differential in `x_pom` |
| `charm_fixed_qp` | the charm quark at fixed `q+`, no photon flux |
| `charm_fixed_qp_xpom` | the same, differential in `x_pom` |
| `scan_flux` | the photon flux alone |
| `saturation_scale` | the saturation scale of a dipole file |

The run scripts loop over `pD0`, `y` and the dipole files. Their settings are environment variables, and the results go to `output/`:

```bash
CHANNEL=0n0n FRAG_TYPE=BCFY ./run_scripts/run_nucleus.sh
Y_VALS="0.0 1.0" PT_VALS="1.0 2.0 4.0" CALLS=1e4 ./run_scripts/run_nucleus.sh
```

## Inputs

Dipoles, MVe parametrization from [rcbkdipole](https://github.com/hejajama/rcbkdipole):

- proton: `data/proton/mve.dat`
- nucleus: `data/Pb/mve/glauber_mve_<b_d>`, one Glauber-sampled file per impact parameter `b_d`

In `input/`:

| Folder | Contents |
|---|---|
| `BCFY_EKO/`, `KK_EKO/` | BCFY and Kniehl-Kramer fragmentation functions, DGLAP-evolved with [eko](https://github.com/NNPDF/eko) in [EKO-FF](https://github.com/patgies/EKO-FF) |
| `HymnD/` | HymnD fragmentation functions, members 0–100 |
| `WS_photon_flux/` | hadronic survival factor `Gamma_AA(b)` for `sigma_NN = 92` mb (default) and `90.85` mb, made with [make_gamma_aa.py](src/make_gamma_aa.py) |
| `Starlight_photon_flux/` | photon flux tables of Petja Paakkinen (EMD probabilities from Starlight v313), for `FLUX_MODEL=STARLIGHT` and `flux_comparison.py` |
| `BK/` | posterior samples of the BK initial condition (see [Uncertainty bands](#uncertainty-bands)) |
| `inclusiveD0/` | inclusive D0 cross section from [inclusive-D0-UPC](https://github.com/patgies/inclusive-D0-UPC), only for `ratio.py` |

`run_scripts/setup_bk_posterior_links.sh` links the BK samples as `data/Pb/bk_posterior/member_NNNN/`. The proton dipoles of the same samples are in `data/proton/bk_posterior/`.

## Settings

| Variable | Values |
|---|---|
| `CHANNEL` | `An0n` (default), `Xn0n`, `0n0n`, `AnAn`, `PL(AnAn)` |
| `FLUX_MODEL` | `EFF` (default) [1]; `STARLIGHT` [1], for `AnAn`, `An0n` and `0n0n`; `PL` (point-like nucleus); `WS` (Woods-Saxon charge distribution) |
| `FRAG_TYPE` | `HymnD` (default) [2], `BCFY` [3], `KniehlKramer` [4] |
| `SCALE_FACTOR` | fragmentation scale `Q = SCALE_FACTOR * mT` (default 1) |

`EFF` is Eq. (4) of [1]. `STARLIGHT` uses the flux tables of the same paper.

`BCFY` and `KniehlKramer` are evolved with DGLAP from their input at the scale `mc`. The grids are in `input/BCFY_EKO/` and `input/KK_EKO/`.

How the photon flux and the fragmentation functions work, and how to change them, is in [notes_photon_flux.pdf](notes_photon_flux.pdf) and [notes_FF.pdf](notes_FF.pdf).

## Uncertainty bands

**Scale uncertainty.** The fragmentation scale `Q` and the renormalization scale of `alpha_s` are varied by factors `0.5` and `2` around `mT`. The band is the envelope of the 7 combinations where the two scales differ by at most a factor of 2.

- `Q`: run `run_scale_variation.sh`, once per fragmentation function.
- `alpha_s`: done in the plotting scripts, and it only changes the diffractive cross section.
- [PbPb_bins.py](plotting_scripts/PbPb_bins.py) and [pPb_bins.py](plotting_scripts/pPb_bins.py) show the band for all three fragmentation functions. [D0.py](plotting_scripts/D0.py) and [D0_sum.py](plotting_scripts/D0_sum.py) show it only for HymnD.

**Other uncertainties**, run with `run_members.sh`:

- `MEMBER_SET=HymnD`: one run per HymnD replica.
- `MEMBER_SET=bk`: one run per posterior sample of the BK fit, with the fragmentation function fixed. The samples are in `input/BK/` and `data/Pb/bk_posterior/`.

When these results are in `output/`, [D0.py](plotting_scripts/D0.py) and [D0_sum.py](plotting_scripts/D0_sum.py) add them in quadrature to the scale uncertainty. The legend says which band is shown. [RpA_fixed_qp.py](plotting_scripts/RpA_fixed_qp.py) shows the BK uncertainty of `R_pA` for the charm quark, from `output/charm_bk/`.

## Normalization

- The nuclear runs give one result per `b_d`. The plotting scripts integrate over `b_d` (Simpson's rule, times `2πb_d`) and multiply by the prefactors.
- For a proton target, the `b_d` integral is replaced by the proton normalization `16.36` mb of MVe.
- The nuclear dipole files stop at `b_d = 32 GeV⁻¹`. The plotting scripts add the dilute edge beyond it: `nuclear_edge_factor` (in [D0.py](plotting_scripts/D0.py)) times the proton-target cross section from `output/proton_baseline/`, made with `run_proton.sh`. `NUCLEAR_EDGE=0` switches it off.
- Units are GeV throughout.

## References

1. K. J. Eskola, V. Guzey, I. Helenius, P. Paakkinen, and H. Paukkunen, "Spatial resolution of dijet photoproduction in near-encounter ultraperipheral nuclear collisions," Phys. Rev. C 110, 054906 (2024) [arXiv:2404.09731].
2. Epele, Hekhorn, Helenius, Paukkunen, and Zurita, "Towards new D meson fragmentation functions" [arXiv:2609.10327].
3. E. Braaten, K.-m. Cheung, S. Fleming, and T.-C. Yuan, "Perturbative QCD fragmentation functions as a model for heavy quark fragmentation," Phys. Rev. D 51, 4819 (1995) [arXiv:hep-ph/9409316].
4. B. A. Kniehl and G. Kramer, "Charmed-hadron fragmentation functions from CERN LEP1 revisited," Phys. Rev. D 74, 037502 (2006) [arXiv:hep-ph/0607306].
