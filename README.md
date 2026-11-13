# Diffractive D0 photoproduction

Code to compute exclusive and diffractive D0 photoproduction `dσ / (dy d²p_D0)` in ultraperipheral collisions in the CGC, for proton and nuclear targets.

If you use it, please cite:

P. Caucal, P. Gimeno-Estivill, E. Iancu, T. Lappi, and F. Salazar, *TMD factorization in diffractive heavy quark production in photon-nucleus collisions*, JHEP 09 (2026) 169 [[arXiv:2606.04169](https://arxiv.org/abs/2606.04169)].

## Build

Needs CMake, GSL and OpenMP.

```bash
mkdir build && cd build
cmake ..
make
```

## Run

```bash
./build/bin/D0 <dipole_file> <pD0> <y> [<FRAG_TYPE>] [<CHANNEL>]
./build/bin/D0 data/proton/mve.dat 3 0 HymnD An0n
```

It prints `pD0`, the exclusive and the diffractive cross section (without prefactors).

Other executables: `D0_xpom` (differential in `x_pom`), `charm_fixed_qp` (charm quark at fixed `q+`, no flux), `charm_fixed_qp_xpom`, `scan_flux` (photon flux only) and `saturation_scale`.

The run scripts are in `run_scripts/`. The settings are environment variables, with the defaults in `run_scripts/config.sh`. Results go to `output/`.

```bash
CHANNEL=0n0n FRAG_TYPE=BCFY ./run_scripts/run_nucleus.sh
Y_VALS="0.0 1.0" PT_VALS="1.0 2.0 4.0" CALLS=1e4 ./run_scripts/run_nucleus.sh
```

The plots are made with the scripts in `plotting_scripts/`.

## Inputs

Dipoles (MVe, from [rcbkdipole](https://github.com/hejajama/rcbkdipole)):

- proton: `data/proton/mve.dat`
- Pb: `data/Pb/mve/glauber_mve_<b_d>`, one file per impact parameter `b_d`

In `input/`:

- `BCFY_EKO/`, `KK_EKO/`: BCFY and Kniehl-Kramer, evolved with [eko](https://github.com/NNPDF/eko) in [EKO-FF](https://github.com/patgies/EKO-FF)
- `HymnD/`: HymnD members 0–100
- `WS_photon_flux/`: `Gamma_AA(b)` tables (made with `src/make_gamma_aa.py`)
- `Starlight_photon_flux/`: flux tables from Petja Paakkinen (Starlight v313)
- `BK/`: BK posterior samples, linked into `data/Pb/bk_posterior/` with `run_scripts/setup_bk_posterior_links.sh`
- `inclusiveD0/`: inclusive results from [inclusive-D0-UPC](https://github.com/patgies/inclusive-D0-UPC), only used in `ratio.py`

## Options

- `CHANNEL`: `An0n` (default), `Xn0n`, `0n0n`, `AnAn`, `PL(AnAn)`
- `FLUX_MODEL`: `EFF` (default, Eq. (4) of [arXiv:2404.09731](https://arxiv.org/abs/2404.09731)), `STARLIGHT`, `PL`, `WS`
- `FRAG_TYPE`: `HymnD` (default, [arXiv:2609.10327](https://arxiv.org/abs/2609.10327)), `BCFY` ([hep-ph/9409316](https://arxiv.org/abs/hep-ph/9409316)), `KniehlKramer` ([hep-ph/0607306](https://arxiv.org/abs/hep-ph/0607306))
- `SCALE_FACTOR`: fragmentation scale `Q = SCALE_FACTOR * mT`

More details on the flux and the fragmentation functions are in [notes_photon_flux.pdf](notes_photon_flux.pdf) and [notes_FF.pdf](notes_FF.pdf).

## Uncertainties

- Scale: `Q` and `mu_R` varied by 0.5 and 2 around `mT` (7-point envelope). `Q` with `run_scale_variation.sh`, `mu_R` in the plotting scripts (only changes the diffractive part).
- HymnD replicas and BK initial condition: `MEMBER_SET=HymnD` or `MEMBER_SET=bk` with `run_members.sh`. `D0.py` and `D0_sum.py` add them in quadrature when they are in `output/`.

## Normalization

The nuclear runs are per `b_d`, and the plotting scripts integrate over `b_d`. For the proton, the `b_d` integral is replaced by `16.36` mb. The Pb dipoles stop at `b_d = 32 GeV⁻¹`; the edge beyond that is added in the plotting scripts using the proton results in `output/proton_baseline/` (`run_proton.sh`). `NUCLEAR_EDGE=0` turns it off.

Units are GeV.
