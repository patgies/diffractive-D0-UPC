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
    "axes.labelsize": 22,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 8,
    "ytick.major.size": 8,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
})

# Large-P tail of exclusive charm photoproduction (Q^2 = 0, z = 1/2) for different dipole amplitudes:
#   F(P) = m^2 [ (z^2 + zbar^2) I_1^2 + I_0^2 ],
#   I_1 = int dr r J1(P r) K1(m r) N(r),  I_0 = int dr r J0(P r) K0(m r) N(r).
# F is the fixed-q+ exclusive cross section without its prefactor (as in Eq. (4.2) of the large-p notes).
# The r integrals are done with Gauss-Legendre between consecutive zeros of the Bessel J (double precision).
# Output: ../plots/exclusive_tail.pdf

M, Z = 1.5, 0.5
QS02, LAM, EC = 0.06, 0.241, 18.9      # MVe fit (data/proton/mve.dat header)
P_PLOT_MAX = 100.0                      # plotted range (the evolved BK curves are noisy beyond it)
SQRT_S = 5360.0                        # as in output/charm/charm_fixed_qp_y00.dat
R_MAX = 40.0                           # K_0, K_1 ~ exp(-m r): negligible beyond
GL_X, GL_W = np.polynomial.legendre.leggauss(24)


def n_mve(r, ec):
    return -np.expm1(-(r**2 * QS02) / 4 * np.log(1 / (r * LAM) + ec * np.e))


def n_gbw(r):
    # GBW with the same Q_s as MVe in the region that matters (log evaluated at r = 1/m)
    return -np.expm1(-(r**2 * QS02) / 4 * np.log(M / LAM + EC * np.e))


def load_bk(filename="../data/proton/mve.dat"):
    """BK solution of AmplitudeLib's format: {y: N(r)} on the log r grid."""
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
    """N(r) at x: linear in the evolution rapidity y = ln(x0/x) between blocks (as AmplitudeLib), then a
    cubic spline in (ln r, ln N)."""
    y = max(np.log(x0 / x), 0.0)
    ys = np.array(sorted(blocks))
    i = min(np.searchsorted(ys, y, side="right") - 1, len(ys) - 2)
    t = (y - ys[i]) / (ys[i + 1] - ys[i])
    n = np.clip((1 - t) * blocks[ys[i]] + t * blocks[ys[i + 1]], 1e-300, None)
    spline = CubicSpline(np.log(r_grid), np.log(n))
    return lambda r: np.exp(spline(np.log(r))), y


def hankel(P, order, amplitude):
    """int_0^R_MAX dr r J_order(P r) K_order(m r) N(r), between consecutive zeros of J_order(P r)."""
    nz = int(R_MAX * P / np.pi) + 2
    zeros = jn_zeros(order, nz) / P if order == 1 else jn_zeros(0, nz) / P
    edges = np.concatenate([[0.0], zeros[zeros < R_MAX], [R_MAX]])
    a, b = edges[:-1, None], edges[1:, None]
    r = 0.5 * (b - a) * GL_X[None, :] + 0.5 * (b + a)
    w = 0.5 * (b - a) * GL_W[None, :]
    J, K = (j1, k1) if order == 1 else (j0, k0)
    return np.sum(w * r * J(P * r) * K(M * r) * amplitude(r))


def spectrum(Ps, amplitude):
    i1 = np.array([hankel(P, 1, amplitude) for P in Ps])
    i0 = np.array([hankel(P, 0, amplitude) for P in Ps])
    return M**2 * ((Z**2 + (1 - Z)**2) * i1**2 + i0**2)


def main():
    Ps = np.geomspace(1.0, 1000.0, 70)
    r_grid, blocks = load_bk()
    bk_x0, _ = bk_amplitude(0.01, r_grid, blocks)
    bk_3, _ = bk_amplitude(1e-3, r_grid, blocks)
    bk_5, _ = bk_amplitude(1e-5, r_grid, blocks)

    curves = [
        (r"GBW", n_gbw, "tab:gray", "-"),
        (r"MV ($e_c=1$)", lambda r: n_mve(r, 1.0), "tab:green", "-"),
        (r"MVe ($e_c=18.9$)", lambda r: n_mve(r, EC), "tab:blue", "-"),
        (r"BK, $x_P=x_0$", bk_x0, "tab:cyan", (0, (4, 3))),
        (r"BK, $x_P=10^{-3}$", bk_3, "tab:orange", "-"),
        (r"BK, $x_P=10^{-5}$", bk_5, "tab:red", "-"),
    ]
    results = {}
    for label, amp, color, ls in curves:
        results[label] = spectrum(Ps, amp)

    fig, (top, bottom) = plt.subplots(2, 1, figsize=(8, 10), sharex=True,
                                      gridspec_kw={"height_ratios": [1.6, 1], "hspace": 0.05})
    keep = Ps <= P_PLOT_MAX * 1.0001
    for label, amp, color, ls in curves:
        F = results[label]
        slope = np.gradient(np.log(F), np.log(Ps))
        lw = 1.6 if label.startswith("BK, $x_P=x_0$") else 2
        top.plot(Ps[keep], F[keep], color=color, ls=ls, lw=lw, label=label)
        bottom.plot(Ps[keep], slope[keep], color=color, ls=ls, lw=lw)
    # The code at fixed q+ = 2p+ (y = 0, proton): x_P = M^2/W^2 grows with P (output/charm, z m^2 g = F/2)
    code = np.loadtxt("../output/charm/charm_fixed_qp_y00.dat", comments="#")
    pc, Fc = code[:, 0], 2 * code[:, 1]
    top.plot(pc, Fc, "ko", ms=5, zorder=5, label=r"BK, $x_P=2m_T/\sqrt{s}$")
    bottom.plot(pc, np.gradient(np.log(Fc), np.log(pc)), "ko", ms=5, zorder=5)
    Fv_at = np.array([spectrum([P], bk_amplitude(2 * np.hypot(P, M) / SQRT_S, r_grid, blocks)[0])[0]
                      for P in pc])  # check: BK at the code's x_P(P) = 2 m_T/sqrt(s) reproduces the dots
    sel = pc >= 2
    print("BK at x_P(P) / code, P >= 2 GeV: min %.5f max %.5f" % ((Fv_at / Fc)[sel].min(), (Fv_at / Fc)[sel].max()))
    print("x_P(P):", " ".join(f"P={P:g}:{2*np.hypot(P, M)/SQRT_S:.2e}" for P in (2, 5, 10, 20, 50, 100)))
    for n in (6, 8, 10):
        bottom.axhline(-n, color="0.6", lw=0.8, ls=":")
        bottom.text(1.08, -n + 0.15, rf"$P_\perp^{{-{n}}}$", fontsize=15, color="0.4")

    top.set_xscale("log")
    top.set_yscale("log")
    top.set_ylim(1e-17, 3e-2)
    top.set_ylabel(r"$F(P_\perp)$ [GeV$^{-2}$]")
    handles, labels = top.get_legend_handles_labels()
    fixed = [i for i, l in enumerate(labels) if not l.startswith("BK, $x_P=2m_T")]
    moving = [i for i, l in enumerate(labels) if i not in fixed]
    leg1 = top.legend([handles[i] for i in fixed], [labels[i] for i in fixed], title=r"fixed $x_P$",
                      fontsize=14, title_fontsize=14, frameon=False, loc="upper right")
    top.add_artist(leg1)
    top.legend([handles[i] for i in moving], [labels[i] for i in moving], title=r"$q^+=2p^+$, $y=0$, proton",
               fontsize=14, title_fontsize=14, frameon=False, loc="lower left")
    bottom.set_ylabel(r"$d\ln F/d\ln P_\perp$")
    bottom.set_xlabel(r"$P_\perp$ [GeV]")
    bottom.set_ylim(-11, -2)
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
