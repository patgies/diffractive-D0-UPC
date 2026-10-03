import argparse
import datetime
import numpy as np

# Gamma_AA(b) = exp[-sigma_NN T_AA(b)], optical Glauber. sigma_NN: 92 mb, 90.8533 mb (Starlight).


HBARC = 0.197327


def thickness(R, a, A, s):
    """T_A(s) [fm^-2]: integral of the WS density along z. int d2s T_A = A."""
    z = np.linspace(-40.0, 40.0, 4001)
    TA = np.array([np.trapezoid(1.0 / (1.0 + np.exp((np.sqrt(z**2 + ss**2) - R) / a)), z) for ss in s])
    return TA * A / np.trapezoid(2 * np.pi * s * TA, s)


def overlap(s, TA, b_fm):
    """T_AA(b) [fm^-2] on a polar grid centred on one nucleus."""
    sp = np.linspace(0.0, 25.0, 1001)
    ph = np.linspace(0.0, np.pi, 721)
    SP, PH = np.meshgrid(sp, ph, indexing="ij")
    TAs = np.interp(sp, s, TA)
    out = []
    for b in b_fm:
        d = np.sqrt(SP**2 + b**2 - 2 * SP * b * np.cos(PH))
        inner = 2 * np.trapezoid(np.interp(d, s, TA, right=0.0), ph, axis=1)
        out.append(np.trapezoid(sp * TAs * inner, sp))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sigma", type=float, required=True, help="sigma_NN [mb]")
    ap.add_argument("--out", required=True)
    ap.add_argument("--R", type=float, default=6.49, help="WS radius [fm]")
    ap.add_argument("--a", type=float, default=0.54, help="WS diffuseness [fm]")
    ap.add_argument("--A", type=float, default=208.0)
    ap.add_argument("--bmax", type=float, default=150.0, help="largest b [GeV^-1]")
    ap.add_argument("--nb", type=int, default=100)
    args = ap.parse_args()

    b = np.linspace(0.0, args.bmax, args.nb)
    s = np.linspace(0.0, 40.0, 4001)
    TAA = overlap(s, thickness(args.R, args.a, args.A, s), b * HBARC)
    gamma = np.exp(-args.sigma * 0.1 * TAA)
    gamma[b == 0.0] = 0.0

    header = "\n".join([
        "Gamma_AA(b) = exp[-sigma_NN * T_AA(b)]  (optical Glauber, Pb+Pb)",
        f"generated      : {datetime.datetime.now():%Y-%m-%d %H:%M:%S} by src/make_gamma_aa.py",
        f"sigma_NN       : {args.sigma:g} mb",
        f"Woods-Saxon    : R = {args.R:g} fm, a = {args.a:g} fm, A = {args.A:g}",
        f"hbar*c         : {HBARC} GeV fm",
        f"b grid         : {args.nb} points, 0 to {args.bmax:g} GeV^-1 (linear)",
        "columns        : b [GeV^-1]   Gamma_AA(b)",
    ])
    np.savetxt(args.out, np.c_[b, gamma], fmt="%.18e", header=header, comments="# ")
    print(f"Saved: {args.out}  (sigma_NN = {args.sigma:g} mb)")


if __name__ == "__main__":
    main()
