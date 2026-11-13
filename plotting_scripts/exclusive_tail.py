import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.special import j0, j1, k0, k1, jn_zeros
from scipy.interpolate import CubicSpline

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 28,
    "axes.titlesize": 26,
    "xtick.labelsize": 25,
    "ytick.labelsize": 25,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 10,
    "ytick.major.size": 10,
    "xtick.minor.size": 5,
    "ytick.minor.size": 5,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.color": "gray",
    "ytick.color": "gray",
    "xtick.labelcolor": "black",
    "ytick.labelcolor": "black",
    "xtick.major.pad": 8,
    "ytick.major.pad": 4,
})

# Large-P tail of exclusive charm photoproduction (Q^2 = 0, z = 1/2) for several dipoles:
#   F(P) = m^2 [ (z^2 + zbar^2) I_1^2 + I_0^2 ],
#   I_1 = int dr r J1(P r) K1(m r) N(r),  I_0 = int dr r J0(P r) K0(m r) N(r).
# F: fixed-q+ exclusive cross section without prefactor. Output: ../plots/exclusive_tail.pdf

M, Z = 1.5, 0.5
QS02, LAM, EC = 0.06, 0.241, 18.9      # MVe fit
P_PLOT_MAX = 100.0                      # plotted range
SQRT_S = 5360.0
R_MAX = 40.0                           # K_0, K_1 negligible beyond
GL_X, GL_W = np.polynomial.legendre.leggauss(24)
# Plotted: alpha_em Nc e_c^2 S_perp/(2 pi^2) * z F [mb/GeV^2], S_perp = sigma0/2
ALPHA_EM, NC, E_C, SIGMA0_MB = 1 / 137, 3, 2 / 3, 16.36
XS = ALPHA_EM * NC * E_C**2 * SIGMA0_MB / (2 * np.pi**2) * Z


def n_mve(r, ec):
    return -np.expm1(-(r**2 * QS02) / 4 * np.log(1 / (r * LAM) + ec * np.e))


def n_gbw(r):
    # GBW with the Q_s of MVe
    return -np.expm1(-(r**2 * QS02) / 4 * np.log(M / LAM + EC * np.e))


def load_bk(filename="../data/proton/mve.dat"):
    """BK solution in AmplitudeLib format: {y: N(r)}."""
    lines = [l for l in open(filename) if not l.startswith("#") or l.startswith("###")]
    rmin, mult, npts = float(lines[0][3:]), float(lines[1][3:]), int(float(lines[2][3:]))
    r = rmin * mult**np.arange(npts)
    blocks, i = {}, 4
    while i + npts < len(lines) + 1 and lines[i].startswith("###"):
        y = float(lines[i][3:])
        blocks[y] = np.array([float(v) for v in lines[i + 1:i + 1 + npts]])
        i += 1 + npts
    return r, blocks


def bk_amplitude(x, r_grid, blocks, x0=0.01):
    """N(r) at x: linear in y = ln(x0/x), spline in (ln r, ln N)."""
    y = max(np.log(x0 / x), 0.0)
    ys = np.array(sorted(blocks))
    i = min(np.searchsorted(ys, y, side="right") - 1, len(ys) - 2)
    t = (y - ys[i]) / (ys[i + 1] - ys[i])
    n = np.clip((1 - t) * blocks[ys[i]] + t * blocks[ys[i + 1]], 1e-300, None)
    spline = CubicSpline(np.log(r_grid), np.log(n))
    return lambda r: np.exp(spline(np.log(r))), y


def hankel(P, order, amplitude):
    """int_0^R_MAX dr r J_order(P r) K_order(m r) N(r)."""
    nz = int(R_MAX * P / np.pi) + 2
    zeros = jn_zeros(order, nz) / P if order == 1 else jn_zeros(0, nz) / P
    edges = np.concatenate([[0.0], zeros[zeros < R_MAX], [R_MAX]])
    a, b = edges[:-1, None], edges[1:, None]
    r = 0.5 * (b - a) * GL_X[None, :] + 0.5 * (b + a)
    w = 0.5 * (b - a) * GL_W[None, :]
    J, K = (j1, k1) if order == 1 else (j0, k0)
    return np.sum(w * r * J(P * r) * K(M * r) * amplitude(r))


def spectrum(Ps, amplitude, parts=False):
    """F(P); parts=True also returns the Q_0 and Q_1 terms."""
    i1 = np.array([hankel(P, 1, amplitude) for P in Ps])
    i0 = np.array([hankel(P, 0, amplitude) for P in Ps])
    q0, q1 = M**2 * i0**2, M**2 * (Z**2 + (1 - Z)**2) * i1**2
    return (q0, q1) if parts else q0 + q1


def main():
    Ps = np.geomspace(1.0, 1000.0, 70)
    r_grid, blocks = load_bk()
    bk_x0, _ = bk_amplitude(0.01, r_grid, blocks)
    bk_3, _ = bk_amplitude(1e-3, r_grid, blocks)
    bk_5, _ = bk_amplitude(1e-5, r_grid, blocks)

    curves = [
        (r"GBW", n_gbw, "0.45", "-"),
        (r"MVe ($e_c=18.9$)", lambda r: n_mve(r, EC), "#2166ac", "-"),
        (r"BK, $x_P=10^{-3}$", bk_3, "#b2182b", "--"),
    ]
    table_only = [
        (r"MV ($e_c=1$)", lambda r: n_mve(r, 1.0)),
        (r"BK, $x_P=x_0$", bk_x0),
        (r"BK, $x_P=10^{-5}$", bk_5),
    ]
    results = {}
    for label, amp, *_ in curves + table_only:
        results[label] = spectrum(Ps, amp)

    fig, (top, bottom) = plt.subplots(2, 1, figsize=(8.5, 11), sharex=True,
                                      gridspec_kw={"height_ratios": [1.6, 1], "hspace": 0.05})
    keep = Ps <= P_PLOT_MAX * 1.0001
    for label, amp, color, ls in curves:
        F = results[label]
        slope = np.gradient(np.log(F), np.log(Ps))
        top.plot(Ps[keep], XS * F[keep], color=color, ls=ls, lw=2.4, label=label)
        bottom.plot(Ps[keep], slope[keep], color=color, ls=ls, lw=2.4)
    # the two terms of the MVe curve
    mve_parts = spectrum(Ps, lambda r: n_mve(r, EC), parts=True)
    for F, color, label in zip(mve_parts, ("#4dac26", "#e08214"), (r"MVe, $Q_0$ term", r"MVe, $Q_1$ term")):
        top.plot(Ps[keep], XS * F[keep], color=color, ls=(0, (1, 0.8)), lw=2, label=label)
        bottom.plot(Ps[keep], np.gradient(np.log(F), np.log(Ps))[keep], color=color, ls=(0, (1, 0.8)), lw=2)
    # Massless quarks at Q^2 > 0: only the Q_1 term remains. GBW massless: P^-10 tail.
    gbw_q1 = spectrum(Ps, n_gbw, parts=True)[1]
    top.plot(Ps[keep], XS * gbw_q1[keep], color="0.45", ls=(0, (1, 0.8)), lw=2, label=r"GBW, $m=0$")
    bottom.plot(Ps[keep], np.gradient(np.log(gbw_q1), np.log(Ps))[keep], color="0.45", ls=(0, (1, 0.8)), lw=2)
    print("MVe Q_0/F:", " ".join(f"P={Ps[i]:.1f}:{(mve_parts[0] / sum(mve_parts))[i]:.2f}"
                                 for i in [np.argmin(abs(Ps - p)) for p in (1, 2, 5, 10, 20, 50, 100)]))
    # The code at fixed q+ = 2p+ (y = 0, proton), output/charm
    code = np.loadtxt("../output/charm/charm_fixed_qp_y00.dat", comments="#")
    pc, Fc = code[:, 0], 2 * code[:, 1]
    top.plot(pc, XS * Fc, "ko", ms=6, zorder=5, label=r"BK, $x_P=2m_T/\sqrt{s}$")
    bottom.plot(pc, np.gradient(np.log(Fc), np.log(pc)), "ko", ms=6, zorder=5)
    Fv_at = np.array([spectrum([P], bk_amplitude(2 * np.hypot(P, M) / SQRT_S, r_grid, blocks)[0])[0]
                      for P in pc])  # check: BK at the code's x_P(P)
    sel = pc >= 2
    print("BK at x_P(P) / code, P >= 2 GeV: min %.5f max %.5f" % ((Fv_at / Fc)[sel].min(), (Fv_at / Fc)[sel].max()))
    print("x_P(P):", " ".join(f"P={P:g}:{2*np.hypot(P, M)/SQRT_S:.2e}" for P in (2, 5, 10, 20, 50, 100)))
    for n in (6, 8, 10):
        bottom.axhline(-n, color="0.6", lw=0.8, ls=":")
        bottom.text(1.08, -n + 0.2, rf"$p_\perp^{{-{n}}}$", fontsize=20, color="0.4")

    top.set_title(r"Exclusive $\gamma+p\rightarrow c\bar c+p$ at $z=0.5$", pad=22, fontsize=28)
    top.set_xscale("log")
    top.set_yscale("log")
    top.set_ylim(1e-19, 1e-4)
    top.set_ylabel(r"$d\sigma/dy\,d^2\mathbf{p}_\perp$ [mb/GeV$^2$]", labelpad=14)
    handles, labels = top.get_legend_handles_labels()
    order = sorted(range(len(labels)), key=lambda i: (0 if labels[i].startswith("GBW") else 1 if labels[i].startswith("MVe")
                                                      else 3 if labels[i].startswith("BK, $x_P=2m_T") else 2))
    top.legend([handles[i] for i in order], [labels[i] for i in order], fontsize=20, frameon=False, loc="lower left")
    bottom.set_ylabel(r"$d\ln\sigma/d\ln p_\perp$", labelpad=14)
    bottom.set_xlabel(r"$p_\perp$ [GeV]", labelpad=12)
    bottom.set_ylim(-11, -2)
    bottom.set_yticks([-10, -8, -6, -4, -2])
    bottom.set_xlim(1, P_PLOT_MAX)
    bottom.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))

    plt.savefig("../plots/exclusive_tail.pdf", bbox_inches="tight")
    print("Saved: ../plots/exclusive_tail.pdf")

    for label in results:
        F = results[label]
        slope = np.gradient(np.log(F), np.log(Ps))
        sel = [np.argmin(abs(Ps - p)) for p in (5, 10, 20, 50, 100, 300, 1000)]
        print(f"{label:20s}", " ".join(f"P={Ps[i]:6.1f}:{slope[i]:6.2f}" for i in sel))


if __name__ == "__main__":
    main()
