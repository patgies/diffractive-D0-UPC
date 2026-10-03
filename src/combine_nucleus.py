import glob
import math
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson
sys.path.insert(0, os.path.dirname(__file__))
from alphas_running import alphas_run

# Integrates the run_nucleus.sh output over the impact parameter b_d
# (Simpson rule, with weight b_d) and multiplies by the prefactor and by 2*pi*pD0.

alphae = 1/137
mc     = 1.5
e_c    = 2/3
Nc     = 3

FMGEV = 5.068
GEVSQR_TO_NB = 1.0e7 / (FMGEV * FMGEV)
GEVSQR_TO_MB = GEVSQR_TO_NB * 1e-6

NUCLEUS = os.environ.get("NUCLEUS", "Pb")
FRAG    = os.environ.get("FRAG_TYPE", "BCFY")
CHANNEL = os.environ.get("CHANNEL", "An0n").translate(str.maketrans('', '', '() '))
# run_nucleus.sh output (flat, or one folder per FRAG) and where the tables go
CENTRAL_DIR = os.environ.get("CENTRAL_DIR", f"../output/{CHANNEL}/central_values")
COMBINED_DIR = os.environ.get("COMBINED_DIR", f"../output/{CHANNEL}/combined")


def read_rapidity(filename):
    with open(filename) as f:
        for line in f:
            if "fixed rapidity y" in line:
                return float(line.split(":")[-1])
    raise ValueError(f"no rapidity header found in {filename}")


def read_data_file(filename):
    """Returns three plain lists (b_d, pD0, dsigma), skipping comment lines."""
    b_d_list, pt_list, dsigma_list = [], [], []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            b_d, pt, dsigma = line.split()
            b_d_list.append(float(b_d))
            pt_list.append(float(pt))
            dsigma_list.append(float(dsigma))
    return b_d_list, pt_list, dsigma_list


def group_by_pt(b_d_list, pt_list, dsigma_list):
    groups = {}
    for b_d, pt, dsigma in zip(b_d_list, pt_list, dsigma_list):
        if pt not in groups:
            groups[pt] = []
        groups[pt].append((b_d, dsigma))
    return groups


def integrate_over_b_d(pairs):
    """Integral of b_d*dsigma(b_d) db_d with Simpson's rule."""
    pairs = sorted(pairs)
    if pairs[0][0] < 0:
        raise ValueError(
            "Negative b_d found -- datafile wasn't a Glauber sample "
            "(e.g. the plain proton dipole has no b_d to integrate over)."
        )
    b_d_values = [pair[0] for pair in pairs]
    weighted = [b_d * dsigma for b_d, dsigma in pairs]
    return simpson(weighted, x=b_d_values)


def prefactor(process, pt):
    """Prefactor of the cross section. No sigma0: the b_d integral already gives the area."""
    if process == "exclusive":
        return alphae * Nc * e_c**2 / (2 * math.pi**2)
    elif process == "diffractive":
        alphas = alphas_run(math.sqrt(pt**2 + mc**2))
        return alphas * alphae * e_c**2 * (Nc**2 - 1) / (8 * math.pi**4)
    raise ValueError(f"unknown process {process}")


def load_results(process):
    name = f"D0_{process}_{FRAG}_{CHANNEL}_{NUCLEUS}_y*.dat"
    results = {}
    for filename in sorted(glob.glob(f"{CENTRAL_DIR}/{name}") or glob.glob(f"{CENTRAL_DIR}/{FRAG}/{name}")):
        y = read_rapidity(filename)
        b_d_list, pt_list, dsigma_list = read_data_file(filename)
        pt_groups = group_by_pt(b_d_list, pt_list, dsigma_list)

        results[y] = []
        for pt, pairs in sorted(pt_groups.items()):
            b_d_integral = integrate_over_b_d(pairs)
            cross_section = (2*math.pi) * b_d_integral * prefactor(process, pt) \
                            * (2*math.pi) * pt * GEVSQR_TO_MB
            results[y].append((pt, cross_section))
    return results


def write_combined(process, results):
    for y, points in sorted(results.items()):
        ytag = str(y).replace(".", "")
        os.makedirs(COMBINED_DIR, exist_ok=True)
        outfile = f"{COMBINED_DIR}/combined_d0_{process}_{FRAG}_{CHANNEL}_{NUCLEUS}_y{ytag}.dat"
        with open(outfile, "w") as f:
            f.write(f"# {process} D0 cross section, {NUCLEUS} target, {FRAG} fragmentation, "
                    f"{CHANNEL} channel, integrated over Glauber b_d\n")
            f.write(f"# fixed rapidity y : {y}\n")
            f.write("# pD0  dsigma_dy_dpD0 [mb/GeV]\n")
            for pt, cs in points:
                f.write(f"{pt:.6g}  {cs:.6g}\n")
        print(f"Wrote {outfile}")


def main():
    results_excl = load_results("exclusive")
    results_diff = load_results("diffractive")
    if not results_excl and not results_diff:
        sys.exit(f"No files found for NUCLEUS={NUCLEUS}, FRAG_TYPE={FRAG}, CHANNEL={CHANNEL} "
                  f"in {CENTRAL_DIR} -- run ../run_scripts/run_nucleus.sh first.")

    write_combined("exclusive", results_excl)
    write_combined("diffractive", results_diff)

    y_values = sorted(set(results_excl) | set(results_diff))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(y_values)))

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for y, c in zip(y_values, colors):
        if y in results_diff:
            points = sorted(results_diff[y])
            pts = [point[0] for point in points]
            cs = [point[1] for point in points]
            ax.plot(pts, cs, color=c, linestyle='-', label=f"diff., $y={y:g}$")
        if y in results_excl:
            points = sorted(results_excl[y])
            pts = [point[0] for point in points]
            cs = [point[1] for point in points]
            ax.plot(pts, cs, color=c, linestyle='--', label=f"excl., $y={y:g}$")

    ax.set_yscale('log')
    ax.set_xlabel(r'$p_{D^0\perp}$ [GeV]')
    ax.set_ylabel(r'$d\sigma/dy\,dp_{D^0\perp}$ (mb/GeV)')
    ax.set_title(f'D0-level $p_{{D^0}}$ spectrum, {NUCLEUS} target ({FRAG}, {CHANNEL})')
    ax.tick_params(axis='both', which='major', labelsize=12)

    handles, labels = ax.get_legend_handles_labels()
    diff_items = [(h, l) for h, l in zip(handles, labels) if l.startswith('diff')]
    excl_items = [(h, l) for h, l in zip(handles, labels) if l.startswith('excl')]
    handles, labels = zip(*(diff_items + excl_items))
    ax.legend(handles, labels, fontsize=7.5, ncol=2)
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    outname = f"../plots/combined_d0_{FRAG}_{CHANNEL}_{NUCLEUS}{NUCLEUS}.png"
    plt.savefig(outname, dpi=300)
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
