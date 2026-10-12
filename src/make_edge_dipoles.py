import argparse
import datetime
import os
import numpy as np

# Dipole files for the dilute edge of the nucleus, glauber_mve_<b_d> with b_d beyond the last
# Glauber sample (b_d = 32 GeV^-1 for data/Pb/mve). There A T_A(b) sigma0/2 < 1 and the nuclear
# amplitude is the proton one times k(b) = A T_A(b) sigma0/2 (arXiv:1309.6963), so each file is the
# proton file with all the amplitudes multiplied by k(b): no BK evolution is needed. With these files
# in the dipole folder, the run scripts and the b_d integral of the plotting scripts include the edge.
#
# usage: python3 src/make_edge_dipoles.py            (Pb: data/Pb/mve/glauber_mve_34 ... _48)
# Remove these files to go back to the Glauber samples only; then the plotting scripts add the edge
# from proton-target runs instead (nuclear_edge_factor in plotting_scripts/D0.py).


HBARC = 0.197327
MB_PER_GEV2 = 0.3894


def thickness(R, a, A, s):
    """A T_A(s) [fm^-2]: integral of the WS density along z. int d2s A T_A = A."""
    z = np.linspace(-40.0, 40.0, 8001)
    TA = np.array([np.trapz(1.0 / (1.0 + np.exp((np.sqrt(z**2 + ss**2) - R) / a)), z) for ss in s])
    return TA * A / np.trapz(2 * np.pi * s * TA, s)


def scaled_file(lines, k, comment):
    """The lines of a BK data file with every amplitude multiplied by k. After the comments (#) come
    four ### lines (minr, multiplier, number of r points, x0) and then, for each rapidity, a ### line
    with y followed by the amplitudes."""
    out = [comment]
    for line in lines:
        if line.startswith("#") or not line.strip():
            out.append(line)
        else:
            out.append(f"{k * float(line):.15e}\n")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proton", default="data/proton/mve.dat", help="proton dipole file")
    parser.add_argument("--outdir", default="data/Pb/mve", help="folder of the nuclear dipole files")
    parser.add_argument("--sigma0", type=float, default=16.36, help="sigma0/2 of the proton fit [mb]")
    parser.add_argument("--A", type=int, default=208)
    parser.add_argument("--R", type=float, default=6.49, help="WS radius [fm] (reproduces the Glauber samples)")
    parser.add_argument("--a", type=float, default=0.542, help="WS diffuseness [fm]")
    parser.add_argument("--b", type=float, nargs="+", default=[34, 36, 38, 40, 42, 44, 46, 48],
                        help="impact parameters b_d [GeV^-1] of the new files")
    args = parser.parse_args()

    with open(args.proton) as f:
        lines = f.readlines()

    s = np.linspace(0.0, 25.0, 5001)
    ATA = thickness(args.R, args.a, args.A, s)
    sigma0 = args.sigma0 / MB_PER_GEV2   # GeV^-2

    for b in args.b:
        k = np.interp(b * HBARC, s, ATA) * HBARC**2 * sigma0   # A T_A(b) sigma0/2
        if k >= 1.0:
            raise SystemExit(f"b_d = {b:g} GeV^-1: A T_A sigma0/2 = {k:.3f} >= 1, not in the dilute edge.")
        comment = (f"# Dilute edge of the nucleus: {args.proton} times k = A T_A(b) sigma0/2 = {k:.6e} at "
                   f"b_d = {b:g} GeV^-1 (A = {args.A}, WS R = {args.R} fm, a = {args.a} fm, sigma0/2 = "
                   f"{args.sigma0} mb), src/make_edge_dipoles.py, {datetime.date.today()}\n")
        outname = os.path.join(args.outdir, f"glauber_mve_{b:g}")
        with open(outname, "w") as f:
            f.writelines(scaled_file(lines, k, comment))
        print(f"{outname}: k = {k:.5f}")


if __name__ == "__main__":
    main()
