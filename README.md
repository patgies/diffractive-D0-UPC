# Diffractive D0 photoproduction

This project computes the exclusive and diffractive D0 photoproduction cross section `dσ / (dy d²p_D0)` in ultraperipheral collisions (UPCs) in the CGC framework.

The code supports different UPC channels, photon fluxes and fragmentation functions, and it can be used for proton and nuclear targets.


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

This prints `pD0`, the exclusive and the diffractive cross section, without the prefactors.

The other executables are `D0_xpom` (differential in `x_pom`), `charm_fixed_qp` (charm quark at fixed `q+`, no photon flux) and `scan_flux` (photon flux alone).

---

## Run scripts

Each calculation has its own script. The results go to `output/`.

| Script | Calculation |
|---|---|
| `run_nucleus.sh` | nuclear target, every `b_d`, `pD0` and `y` |
| `run_scale_variation.sh` | the same with the fragmentation scale `Q = 0.5 mT` and `2 mT` |
| `run_proton.sh` | `TARGET=pA` for p+Pb at 8.16 TeV |
| `run_xpom.sh` | cross section differential in `x_pom` |
| `run_members.sh` | HymnD replicas (`MEMBER_SET=HymnD`) or BK posterior samples (`MEMBER_SET=bk`) |
| `run_charm_fixed_qp.sh` | charm quark at fixed `q+` |
| `run_flux_scan.sh` | photon flux alone |

The settings are environment variables. The defaults are in [run_scripts/config.sh](run_scripts/config.sh).

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

---

## UPC channel and photon flux

The channel is selected through `CHANNEL`:

- `An0n`
- `Xn0n`
- `0n0n`
- `AnAn`
- `PL(AnAn)`

The photon flux is selected through `FLUX_MODEL`: `EFF` (default), `STARLIGHT`, `PL` or `WS`.

The photon flux, and how to use your own, is explained in [notes_photon_flux.pdf](notes_photon_flux.pdf).

---

## Fragmentation functions

The fragmentation function is selected through `FRAG_TYPE`:

- `BCFY`
- `KniehlKramer`
- `HymnD` (default)

The fragmentation scale is `Q = SCALE_FACTOR * mT`.

The fragmentation functions, and how to add your own, are explained in [notes_FF.pdf](notes_FF.pdf).

The HymnD members 0–100 are in `input/HymnD/`.

---

## Differential cross section and normalization

The output of the nuclear runs is per impact parameter `b_d`. The plotting scripts integrate over `b_d` (Simpson's rule, multiplied by `2πb_d`) and multiply by the prefactors.

For a proton target, the impact-parameter integral is replaced by the proton normalization `16.36` mb of the MVe dipole parametrization.

Units are GeV throughout.
