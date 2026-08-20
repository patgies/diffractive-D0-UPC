import glob
import math
import os
import statistics
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.integrate import simpson
sys.path.insert(0, os.path.dirname(__file__))
from alphas_running import alphas_run

# make the plot look nicer (same style as xpom_spectrum.py / fixedW_spectrum.py)
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 20,
    "axes.titlesize": 20,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 8,
    "ytick.major.size": 8,
    "xtick.minor.size": 4,
    "ytick.minor.size": 4,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
})

# This script reads the data files written by ../run_many_nucleus.sh
# (produced by the D0 program), integrates them over the target-nucleus
# impact parameter b using Simpson's rule, and makes one plot with every
# combination of process (exclusive/diffractive) and fragmentation function
# (BCFY/KniehlKramer) on it. Each rapidity gets its own color, and each
# (process, fragmentation) combination gets its own linestyle.

alphae = 1/137
mc     = 1.5          # charm mass in GeV
e_c    = 2/3
Nc     = 3
sigma0 = 16.36         # mb

# conversion factor from GeV^-2 to mb
FMGEV = 5.068
GEVSQR_TO_NB = 1.0e7 / (FMGEV * FMGEV)
GEVSQR_TO_MB = GEVSQR_TO_NB * 1e-6

NUCLEUS = os.environ.get("NUCLEUS", "Pb")
CHANNEL = os.environ.get("CHANNEL", "An0n").translate(str.maketrans('', '', '() '))

PROCESSES = ["diffractive", "exclusive"]
FRAG_TYPES = ["LHAPDF"]
MAX_Y_TO_PLOT = 3.0   # keeps the plot from getting too crowded

# which linestyle to use for each (process, fragmentation) combination
LINESTYLES = {
    ("diffractive", "BCFY"):          "-",
    ("diffractive", "KniehlKramer"):  "--",
    ("diffractive", "LHAPDF"):        "-",
    ("exclusive",   "BCFY"):          ":",
    ("exclusive",   "KniehlKramer"):  "-.",
    ("exclusive",   "LHAPDF"):        ":",
}


def read_rapidity(filename):
    with open(filename) as f:
        for line in f:
            if "fixed rapidity y" in line:
                return float(line.split(":")[-1])
    raise ValueError(f"no rapidity header found in {filename}")


def read_data_file(filename):
    """Returns three plain lists (b, pD0, dsigma), skipping comment lines."""
    b_list = []
    pt_list = []
    dsigma_list = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            b, pt, dsigma = line.split()
            b_list.append(float(b))
            pt_list.append(float(pt))
            dsigma_list.append(float(dsigma))
    return b_list, pt_list, dsigma_list


def group_by_pt(b_list, pt_list, dsigma_list):
    groups = {}
    for b, pt, dsigma in zip(b_list, pt_list, dsigma_list):
        if pt not in groups:
            groups[pt] = []
        groups[pt].append((b, dsigma))
    return groups


def integrate_over_b(pairs):
    """int b*dsigma(b) db via Simpson's rule (the radial Glauber-b measure)."""
    pairs = sorted(pairs)   # sorts by b first since these are (b, dsigma) tuples
    if pairs[0][0] < 0:
        raise ValueError(
            "Negative b found -- datafile wasn't a Glauber sample "
            "(e.g. the plain proton dipole has no b to integrate over)."
        )
    b_values = [pair[0] for pair in pairs]
    weighted = [b * dsigma for b, dsigma in pairs]
    return simpson(weighted, x=b_values)


def prefactor(process, pt):
    """Physical prefactor, matching plot_pt_spectrum.py's convention."""
    if process == "exclusive":
        return alphae * Nc * e_c**2 * sigma0 / (2 * math.pi**2)
    elif process == "diffractive":
        alphas = alphas_run(math.sqrt(pt**2 + mc**2))
        return alphas * alphae * e_c**2 * (Nc**2 - 1) * sigma0 / (8 * math.pi**4)
    raise ValueError(f"unknown process {process}")


def load_results(process, frag, base_dir=".."):
    """Read all files for one (process, frag) combination and return {y: [(pt, cross_section), ...]}."""
    pattern = f"{base_dir}/files/D0_{process}_{frag}_{CHANNEL}_{NUCLEUS}_y*.dat"
    results = {}
    for filename in sorted(glob.glob(pattern)):
        y = read_rapidity(filename)
        b_list, pt_list, dsigma_list = read_data_file(filename)
        pt_groups = group_by_pt(b_list, pt_list, dsigma_list)

        results[y] = []
        for pt, pairs in sorted(pt_groups.items()):
            b_integral = integrate_over_b(pairs)
            # 2*pi*b_integral is the Glauber transverse-plane (b) integral.
            # 2*pi*pt is the Jacobian from d^2pD0 to dpD0.
            cross_section = (2*math.pi) * b_integral * prefactor(process, pt) \
                            * (2*math.pi) * pt * GEVSQR_TO_MB
            results[y].append((pt, cross_section))
    return results


def _band_from_members(process, frag, member_dir_pattern):
    """Shared aggregation for load_lhapdf_band/load_bk_band: combine one
    run_many_nucleus.sh output set per member directory into a mean +/-
    standard-deviation band per rapidity (both uncertainty sources here are
    sampled sets -- LHAPDF replicas, BK posterior samples -- not Hessian
    eigenvectors, so plain standard deviation across members is the right
    prescription, not a Hessian formula).

    Returns {y: (pt_values, means, stds)}, or {} if no member directories
    are found yet.
    """
    member_dirs = sorted(glob.glob(member_dir_pattern))
    if not member_dirs:
        return {}

    per_member_results = [load_results(process, frag, base_dir=member_dir) for member_dir in member_dirs]
    per_member_results = [results for results in per_member_results if results]
    if not per_member_results:
        return {}

    band = {}
    for y in per_member_results[0]:
        pt_values = sorted(pt for pt, _ in per_member_results[0][y])
        cross_sections_by_pt = []
        for pt in pt_values:
            values = [dict(results[y])[pt] for results in per_member_results if y in results]
            cross_sections_by_pt.append(values)

        means = [statistics.mean(values) for values in cross_sections_by_pt]
        stds = [statistics.pstdev(values) for values in cross_sections_by_pt]
        band[y] = (pt_values, means, stds)

    return band


def load_lhapdf_band(process, frag="LHAPDF"):
    """LHAPDF-fragmentation-replica band (see run_lhapdf_members_roihu.sbatch):
    dipole amplitude held fixed at the central data/Pb/mve/ set, fragmentation
    function varied across the 101 replica members. {} until
    run_lhapdf_members_roihu.sbatch's output has been copied back into
    ../files/lhapdf/member_*.
    """
    return _band_from_members(process, frag, "../files/lhapdf/member_*")


def load_bk_band(process, frag="LHAPDF"):
    """BK-initial-condition posterior band (see
    run_bk_posterior_members_roihu.sbatch): fragmentation function held fixed
    at the central LHAPDF member, dipole amplitude varied across the 100
    Bayesian BK-IC posterior samples in bk/. {} until
    run_bk_posterior_members_roihu.sbatch's output has been copied back into
    ../files/bk_posterior/member_*.
    """
    return _band_from_members(process, frag, "../files/bk_posterior/member_*")


def load_scale_band(process, frag="LHAPDF"):
    """LHAPDF fragmentation-scale variation envelope (see
    run_lhapdf_scale_variation.sh): central member only, Q = SCALE_FACTOR *
    mt0 varied at the conventional 0.5x/2x around the central scale
    (SCALE_FACTOR=1.0, i.e. the plain ../files/D0_..._LHAPDF_... run). This
    is a scale-convention envelope (min/max), not a statistical replica/
    posterior sample, so it's built differently from
    load_lhapdf_band/load_bk_band -- min and max across the 3 runs, not
    mean +/- std.

    Returns {y: (pt_values, lower, upper)}, or {} until
    run_lhapdf_scale_variation.sh's output exists in ../files/lhapdf_scale/.
    """
    scale_dirs = sorted(glob.glob("../files/lhapdf_scale/factor_*"))
    if not scale_dirs:
        return {}

    central = load_results(process, frag)
    variations = [load_results(process, frag, base_dir=d) for d in scale_dirs]
    variations = [results for results in variations if results]
    all_results = [central] + variations
    if not central or not variations:
        return {}

    envelope = {}
    for y in central:
        pt_values = sorted(pt for pt, _ in central[y])
        lower = []
        upper = []
        for pt in pt_values:
            values = [dict(results[y])[pt] for results in all_results if y in results and pt in dict(results[y])]
            lower.append(min(values))
            upper.append(max(values))
        envelope[y] = (pt_values, lower, upper)

    return envelope


def main():
    all_results = {}
    for process in PROCESSES:
        for frag in FRAG_TYPES:
            results = load_results(process, frag)
            if results:
                all_results[(process, frag)] = results

    # If run_lhapdf_members_roihu.sbatch's / run_bk_posterior_members_roihu.sbatch's
    # output has been copied back into ../files/lhapdf/ and ../files/bk_posterior/,
    # plot those error bands (each an independent uncertainty source, not
    # combined) alongside/instead of the single-member LHAPDF line.
    lhapdf_bands = {process: load_lhapdf_band(process) for process in PROCESSES}
    bk_bands = {process: load_bk_band(process) for process in PROCESSES}
    scale_bands = {process: load_scale_band(process) for process in PROCESSES}

    if (not all_results and not any(lhapdf_bands.values())
            and not any(bk_bands.values()) and not any(scale_bands.values())):
        sys.exit(f"No files found for NUCLEUS={NUCLEUS}, CHANNEL={CHANNEL} "
                  "in ../files/ -- run ../run_many_nucleus.sh first.")

    # collect every rapidity value that shows up in any of the results,
    # keeping only y <= MAX_Y_TO_PLOT so the plot doesn't get too crowded
    y_values = set()
    for results in all_results.values():
        for y in results:
            if y <= MAX_Y_TO_PLOT:
                y_values.add(y)
    for bands in (lhapdf_bands, bk_bands, scale_bands):
        for band in bands.values():
            for y in band:
                if y <= MAX_Y_TO_PLOT:
                    y_values.add(y)
    y_values = sorted(y_values)

    # ColorBrewer's "Blues" sequential scale (rapidity is ordered: light ->
    # dark = increasing y). A plain hex interpolation within one hue reads as
    # too similar step-to-step; this one is built for perceptual spacing
    # between adjacent steps, not just even numeric spacing.
    blue_ramp = [
        "#deebf7", "#c6dbef", "#9ecae1", "#6baed6",
        "#4292c6", "#2171b5", "#08519c", "#08306b",
    ]
    colors = {}
    for i, y in enumerate(y_values):
        step = round(i * (len(blue_ramp) - 1) / max(len(y_values) - 1, 1))
        colors[y] = blue_ramp[step]

    plt.figure(figsize=(7.5, 6.5))

    for process, frag in LINESTYLES:
        linestyle = LINESTYLES[(process, frag)]

        have_any_band = frag == "LHAPDF" and (lhapdf_bands.get(process) or bk_bands.get(process)
                                               or scale_bands.get(process))
        if have_any_band:
            # LHAPDF-replica band: filled, no hatch. BK-posterior band:
            # hatched outline only. Scale-variation envelope: dashed outline
            # only. Three independent uncertainty sources, kept visually
            # distinguishable when more than one is present at once.
            for band, is_envelope, fill_kwargs in (
                (lhapdf_bands.get(process), False, dict(alpha=0.25, linewidth=0)),
                (bk_bands.get(process), False, dict(alpha=0.0, edgecolor="none", hatch="////", linewidth=0)),
                (scale_bands.get(process), True, dict(alpha=0.0, linestyle="--", linewidth=1.2)),
            ):
                if not band:
                    continue
                for y in sorted(band):
                    if y > MAX_Y_TO_PLOT:
                        continue
                    if is_envelope:
                        pt_values, lower, upper = band[y]
                    else:
                        pt_values, means, stds = band[y]
                        lower = [max(mean - std, 1e-30) for mean, std in zip(means, stds)]
                        upper = [mean + std for mean, std in zip(means, stds)]
                    edge_kwargs = dict(fill_kwargs)
                    edgecolor = edge_kwargs.pop("edgecolor", colors[y])
                    plt.fill_between(pt_values, lower, upper, color=colors[y], edgecolor=edgecolor, **edge_kwargs)
            # Central line: prefer a replica/posterior band's mean (both
            # centered on member-0000); fall back to the plain central run
            # if only the scale envelope is available.
            central_band = lhapdf_bands.get(process) or bk_bands.get(process)
            if central_band:
                for y in sorted(central_band):
                    if y > MAX_Y_TO_PLOT:
                        continue
                    pt_values, means, _ = central_band[y]
                    plt.plot(pt_values, means, color=colors[y], linestyle=linestyle)
            else:
                results = all_results.get((process, frag))
                if results:
                    for y in sorted(results):
                        if y > MAX_Y_TO_PLOT:
                            continue
                        points = sorted(results[y])
                        plt.plot([p[0] for p in points], [p[1] for p in points],
                                 color=colors[y], linestyle=linestyle)
            continue

        results = all_results.get((process, frag))
        if not results:
            continue
        for y in sorted(results):
            if y > MAX_Y_TO_PLOT:
                continue
            points = sorted(results[y])
            pt_values = [pair[0] for pair in points]
            cross_section_values = [pair[1] for pair in points]
            plt.plot(pt_values, cross_section_values, color=colors[y], linestyle=linestyle)

    plt.yscale("log")
    plt.xlabel(r"$p_{D^0}$ [GeV]")
    plt.ylabel(r"$d\sigma/dy\,dp_{D^0}$ [mb/GeV]", labelpad=15)
    plt.title(f"$D^0$ photoproduction, {NUCLEUS}+{NUCLEUS} UPC ({CHANNEL})", pad=15)

    y_handles = [Line2D([0], [0], color=colors[y], linestyle="-", label=f"y={y:g}") for y in y_values]
    y_legend = plt.legend(handles=y_handles, loc="upper right", fontsize=17)
    plt.gca().add_artist(y_legend)

    style_handles = []
    for (process, frag), linestyle in LINESTYLES.items():
        has_lhapdf_band = frag == "LHAPDF" and lhapdf_bands.get(process)
        has_bk_band = frag == "LHAPDF" and bk_bands.get(process)
        has_scale_band = frag == "LHAPDF" and scale_bands.get(process)
        if has_lhapdf_band or has_bk_band or has_scale_band:
            style_handles.append(Line2D([0], [0], color="black", linestyle=linestyle, label=f"{process}, {frag}"))
            if has_lhapdf_band:
                style_handles.append(Patch(facecolor="black", alpha=0.25, label="  ± LHAPDF replica std"))
            if has_bk_band:
                style_handles.append(Patch(facecolor="none", edgecolor="black", hatch="////",
                                            label="  ± BK-IC posterior std"))
            if has_scale_band:
                style_handles.append(Line2D([0], [0], color="black", linestyle="--", linewidth=1.2,
                                             label="  scale var. (0.5x/2x mt0)"))
        elif (process, frag) in all_results:
            style_handles.append(Line2D([0], [0], color="black", linestyle=linestyle, label=f"{process}, {frag}"))
    plt.legend(handles=style_handles, loc="lower left", fontsize=15)

    plt.tight_layout()
    outname = f"../plots/cross_section_{CHANNEL}_{NUCLEUS}.pdf"
    plt.savefig(outname)
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
