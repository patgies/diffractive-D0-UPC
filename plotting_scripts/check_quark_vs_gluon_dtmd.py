import argparse
import math
import numpy as np
from scipy.special import jv, kv
from scipy.interpolate import interp1d

# Check of the quark/gluon DTMD ratio of the diffractive_TMD cross section (Fig. xpom_quark_vs_gluon),
# independent of src/integrands.cpp: the K-integrated gluon and quark DTMD terms are computed directly
# from the formulas of the manuscript (hard factors H and H', gluon DTMD with G, quark DTMD with Q_n)
# at fixed q+ = 2p+ (z = 1/2), for a proton target, at several x = M^2/(x_po W^2) and charm masses.
# The common prefactor alpha_em alpha_s e_q^2 S_perp/(8 pi^4) cancels in the ratio.
# Both terms use x = M_qq^2/(x_po W^2) (M_qg^2 differs by ~2% at p = 8 GeV), and only the DTMDs depend
# on the mass in the comparison (H'/H ~ 3.3 for all masses).
#
# usage: python3 check_quark_vs_gluon_dtmd.py [--p 8] [--y 0] [--masses 1.5 0.3 0.05]

PROTON = "../data/proton/mve.dat"
NC, CF = 3.0, 4.0 / 3.0
Z = 0.5
SS = 5360.0   # sqrt(s_NN) [GeV]
X_VALUES = (0.9, 0.5, 0.2, 0.05)


def read_bk(filename):
    """r grid, rapidities Y and amplitudes N[Y, r] of a BK data file. After the comments (#) come four ###
    lines (minr, multiplier, number of r points, x0) and then, for each Y, a ### line followed by N."""
    head, blocks = [], []
    with open(filename) as f:
        for line in f:
            if line.startswith("###"):
                if len(head) < 4:
                    head.append(float(line[3:]))
                else:
                    blocks.append((float(line[3:]), []))
            elif line.startswith("#") or not line.strip():
                continue
            else:
                blocks[-1][1].append(float(line))
    minr, mult, nr, x0 = head
    r = minr * mult ** np.arange(int(nr))
    return r, np.array([b[0] for b in blocks]), np.array([b[1] for b in blocks]), x0


R_TAB, Y_TAB, N_TAB, X0 = read_bk(PROTON)


def dipole(x_po):
    """N(r) at Y = ln(x0/x_po), frozen at the initial condition above x0 (as XPO_FROZEN in the code)."""
    Y = max(0.0, math.log(X0 / x_po))
    i = min(np.searchsorted(Y_TAB, Y), len(Y_TAB) - 1)
    f = interp1d(np.log(R_TAB), N_TAB[i], bounds_error=False, fill_value=(0.0, 1.0))
    return lambda r: f(np.log(np.maximum(r, 1e-12)))


def radial(n, k, omega, N):
    """omega^2 int dR R J_n(k R) K_n(omega R) N(R)"""
    rmax = min(100.0, 40.0 / omega)   # K_n(omega R) is negligible beyond
    R = np.linspace(1e-6, rmax, 60000)
    return omega**2 * np.trapz(R * jv(n, k * R) * kv(n, omega * R) * N(R), R)


def gluon_dtmd(x, N, kmax):
    """int d2K (Nc^2 - 1)/(1 - x) G^2, with Omega^2 = x K^2/(1 - x) and N_g = 2N - N^2"""
    Ng = lambda R: 2 * N(R) - N(R) ** 2
    K = np.geomspace(0.01, kmax, 160)
    vals = [(NC**2 - 1) / (1 - x) * radial(2, k, math.sqrt(x / (1 - x)) * k, Ng) ** 2 for k in K]
    return np.trapz(2 * np.pi * K * np.array(vals), K)


def quark_dtmd(x, m, N, kmax):
    """int d2K Nc x/(1-x) (K^2+m^2)/(xK^2+m^2) [Q_1^2 + (1-x) m^2/(xK^2+m^2) Q_0^2],
    with Omega^2 = (m^2 + x K^2)/(1 - x) and the quark dipole N"""
    K = np.geomspace(0.01, kmax, 160)
    vals = []
    for k in K:
        omega = math.sqrt((m**2 + x * k**2) / (1 - x))
        q1, q0 = radial(1, k, omega, N), radial(0, k, omega, N)
        xk2m2 = x * k**2 + m**2
        vals.append(NC * x / (1 - x) * (k**2 + m**2) / xk2m2 * (q1**2 + (1 - x) * m**2 / xk2m2 * q0**2))
    return np.trapz(2 * np.pi * K * np.array(vals), K)


def hard_factors(p, m):
    """H (gluon term) and H' (quark term), without alpha_em alpha_s e_q^2"""
    zbar = 1 - Z
    p2, m2 = p**2, m**2
    H = Z * ((Z**2 + zbar**2) * (p2**2 + m2**2) + 2 * m2 * p2) / (p2 + m2) ** 4
    H_quark = 2 * CF * Z * Z * (1 + zbar**2) / p2**2
    return H, H_quark


def main():
    parser = argparse.ArgumentParser(description="Quark/gluon DTMD ratio at fixed q+ = 2p+")
    parser.add_argument("--p", type=float, default=8.0, help="charm transverse momentum [GeV]")
    parser.add_argument("--y", type=float, default=0.0, help="rapidity")
    parser.add_argument("--masses", type=float, nargs="+", default=[1.5, 0.3, 0.05], help="charm masses [GeV]")
    args = parser.parse_args()
    p = args.p

    for m in args.masses:
        mt = math.sqrt(p**2 + m**2)
        w2 = 2 * mt * math.exp(args.y) * SS   # sqrt(2) ss q+ with q+ = 2p+ = sqrt(2) m_T e^y
        Mqq2 = (p**2 + m**2) / (Z * (1 - Z))
        H, H_quark = hard_factors(p, m)
        print(f"# m = {m:g} GeV, p = {p:g} GeV, y = {args.y:g}: H'/H = {H_quark / H:.2f}")
        print("#    x      x_po     gluon      quark    q/g DTMD  q/g with H, H'")
        for x in X_VALUES:
            x_po = Mqq2 / (x * w2)
            N = dipole(x_po)
            g = gluon_dtmd(x, N, 2 * mt)   # K cuts of the manuscript: 2 m_T (gluon), 2 p (quark)
            q = quark_dtmd(x, m, N, 2 * p)
            print(f"  {x:5.2f}  {x_po:.2e}  {g:.3e}  {q:.3e}  {q / g:8.3f}  {q * H_quark / (g * H):8.3f}")


if __name__ == "__main__":
    main()
