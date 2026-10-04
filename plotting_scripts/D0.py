import glob
import math
import os
import statistics
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import ticker
from scipy.integrate import simpson
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))   # alphas_running.py
from alphas_running import alphas_run

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 28,
    "axes.titlesize": 20,
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
    "xtick.color": "black",
    "ytick.color": "black",
    "xtick.labelcolor": "black",
    "ytick.labelcolor": "black",
    "xtick.major.pad": 8,
    "ytick.major.pad": 4,
})

# Reads the D0 output of run_nucleus.sh, integrates it over the impact parameter b_d
# and plots exclusive and diffractive spectra.

alphae = 1/137
mc     = 1.5
e_c    = 2/3
Nc     = 3

FMGEV = 5.068
GEVSQR_TO_NB = 1.0e7 / (FMGEV * FMGEV)
GEVSQR_TO_MB = GEVSQR_TO_NB * 1e-6

NUCLEUS = os.environ.get("NUCLEUS", "Pb")
CHANNEL = os.environ.get("CHANNEL", "An0n").translate(str.maketrans('', '', '() '))

CENTRAL_DIR = os.environ.get("CENTRAL_DIR", f"../output/{CHANNEL}/central_values")
SCALE_DIR = os.environ.get("SCALE_DIR", f"../output/{CHANNEL}/scale_variation/HymnD")

PROCESSES = ["diffractive", "exclusive"]
FRAG_TYPES = ["HymnD"]
MIN_Y_TO_PLOT = -1.0
MAX_Y_TO_PLOT = 2.0

SCALE_FACTOR_TAGS = {0.5: "0.5", 1.0: None, 2.0: "2.0"}
SCALE_FACTORS = list(SCALE_FACTOR_TAGS)
EXPECTED_SCALE_FACTORS = [tag for tag in SCALE_FACTOR_TAGS.values() if tag]


def _mu_f_dir(factor):
    """Folder with the D0_*.dat files for a given mu_F factor (see SCALE_FACTOR_TAGS)."""
    tag = SCALE_FACTOR_TAGS[factor]
    return CENTRAL_DIR if tag is None else f"{SCALE_DIR}/factor_{tag}"

LINESTYLES = {
    ("diffractive", "BCFY"):          "-",
    ("diffractive", "KniehlKramer"):  "--",
    ("diffractive", "HymnD"):        "-",
    ("exclusive",   "BCFY"):          ":",
    ("exclusive",   "KniehlKramer"):  "-.",
    ("exclusive",   "HymnD"):        ":",
}
LINEWIDTH = 2.0


def read_rapidity(filename):
    with open(filename) as f:
        for line in f:
            if "fixed rapidity y" in line:
                return float(line.split(":")[-1])
    raise ValueError(f"no rapidity header found in {filename}")


def read_data_file(filename):
    """Returns three plain lists (b_d, pD0, dsigma), skipping comment lines."""
    b_d_list = []
    pt_list = []
    dsigma_list = []
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


def prefactor(process, pt, mu_r_factor=1.0):
    """Prefactor of the cross section. No sigma0: the b_d integral already gives the area.
    mu_r_factor multiplies the scale of alpha_s. Only the diffractive process uses it."""
    if process == "exclusive":
        return alphae * Nc * e_c**2 / (2 * math.pi**2)
    elif process == "diffractive":
        alphas = alphas_run(mu_r_factor * math.sqrt(pt**2 + mc**2))
        return alphas * alphae * e_c**2 * (Nc**2 - 1) / (8 * math.pi**4)
    raise ValueError(f"unknown process {process}")


def load_results(process, frag, data_dir=None, mu_r_factor=1.0):
    """Read all files for one (process, frag) in data_dir or data_dir/<frag>/
    and return {y: [(pt, cross_section), ...]}."""
    data_dir = CENTRAL_DIR if data_dir is None else data_dir
    name = f"D0_{process}_{frag}_{CHANNEL}_{NUCLEUS}_y*.dat"
    filenames = glob.glob(f"{data_dir}/{name}") or glob.glob(f"{data_dir}/{frag}/{name}")
    results = {}
    for filename in sorted(filenames):
        y = read_rapidity(filename)
        b_d_list, pt_list, dsigma_list = read_data_file(filename)
        pt_groups = group_by_pt(b_d_list, pt_list, dsigma_list)

        MIN_B_D_SAMPLES = 10
        full_b_d_count = max((len(pairs) for pairs in pt_groups.values()), default=0)
        pt_groups = {
            pt: pairs for pt, pairs in pt_groups.items()
            if len(pairs) == full_b_d_count and len(pairs) >= MIN_B_D_SAMPLES
        }

        if not pt_groups:
            continue

        points = []
        for pt, pairs in sorted(pt_groups.items()):
            b_d_integral = integrate_over_b_d(pairs)
            cross_section = (2*math.pi) * b_d_integral * prefactor(process, pt, mu_r_factor) \
                            * (2*math.pi) * pt * GEVSQR_TO_MB
            points.append((pt, cross_section))
        results[y] = points
    return results


def _band_from_members(process, frag, member_dir_pattern):
    """Band from the mean and standard deviation over the member folders.
    Returns {y: (pt_values, means, stds)}, or {} if none are found."""
    member_dirs = sorted(glob.glob(member_dir_pattern))
    if not member_dirs:
        return {}

    per_member_results = [load_results(process, frag, data_dir=member_dir) for member_dir in member_dirs]
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


def load_hymnd_replica_band(process, frag="HymnD"):
    """Band from the HymnD replicas (one fragmentation function per member). {} if there is no data."""
    return _band_from_members(process, frag, f"../output/{CHANNEL}/HymnD_band/member_*")


def load_bk_band(process, frag="HymnD"):
    """Band from the BK posterior samples (one dipole per sample). {} if there is no data."""
    return _band_from_members(process, frag, f"../output/{CHANNEL}/bk_band/member_*")


def _summed_results(data_dir, frag="HymnD", mu_r_factor=1.0):
    """diffractive+exclusive summed per (y, pt) for one directory of D0_*.dat
    files. Returns {y: [(pt, total), ...]}, or {} if neither process has data there.
    """
    diff = load_results("diffractive", frag, data_dir=data_dir, mu_r_factor=mu_r_factor)
    excl = load_results("exclusive", frag, data_dir=data_dir, mu_r_factor=mu_r_factor)
    if not diff and not excl:
        return {}
    totals = {}
    for y in set(diff) | set(excl):
        diff_points = dict(diff.get(y, []))
        excl_points = dict(excl.get(y, []))
        pt_values = sorted(set(diff_points) | set(excl_points))
        totals[y] = [(pt, diff_points.get(pt, 0.0) + excl_points.get(pt, 0.0)) for pt in pt_values]
    return totals


def _summed_band_from_members(member_dir_pattern, frag="HymnD"):
    """Like _band_from_members, but first adds diffractive+exclusive for each member.
    Returns {y: (pt_values, means, stds)}, or {} if none are found."""
    member_dirs = sorted(glob.glob(member_dir_pattern))
    if not member_dirs:
        return {}

    per_member_totals = [_summed_results(member_dir, frag) for member_dir in member_dirs]
    per_member_totals = [totals for totals in per_member_totals if totals]
    if not per_member_totals:
        return {}

    band = {}
    for y in per_member_totals[0]:
        pt_values = sorted(pt for pt, _ in per_member_totals[0][y])
        cross_sections_by_pt = []
        for pt in pt_values:
            values = [dict(totals[y])[pt] for totals in per_member_totals if y in totals]
            cross_sections_by_pt.append(values)

        means = [statistics.mean(values) for values in cross_sections_by_pt]
        stds = [statistics.pstdev(values) for values in cross_sections_by_pt]
        band[y] = (pt_values, means, stds)

    return band


def load_hymnd_replica_sum_band(frag="HymnD"):
    """Band from the HymnD replicas for the diffractive+exclusive sum. {} if there is no data."""
    return _summed_band_from_members(f"../output/{CHANNEL}/HymnD_band/member_*", frag)


def load_bk_sum_band(frag="HymnD"):
    """Band from the BK posterior samples for the diffractive+exclusive sum. {} if there is no data."""
    return _summed_band_from_members(f"../output/{CHANNEL}/bk_band/member_*", frag)


def _scale_combos():
    """The 7-point (mu_F, mu_R) grid: SCALE_FACTORS for both, without the two
    points where mu_F/mu_R is 4 or 1/4 (e.g. arXiv:2506.09893)."""
    return [(mu_f, mu_r) for mu_f in SCALE_FACTORS for mu_r in SCALE_FACTORS
            if 0.5 <= mu_f / mu_r <= 2.0]


def _envelope_from_results(results_by_combo):
    """Minimum and maximum over {(mu_f, mu_r): {y: [(pt, val), ...]}}.
    Returns {y: (pt_values, lower, upper)}."""
    central = results_by_combo[(1.0, 1.0)]
    envelope = {}
    for y in central:
        pt_values = sorted(pt for pt, _ in central[y])
        lower = []
        upper = []
        for pt in pt_values:
            values = [dict(results[y])[pt] for results in results_by_combo.values()
                      if y in results and pt in dict(results[y])]
            lower.append(min(values))
            upper.append(max(values))
        envelope[y] = (pt_values, lower, upper)
    return envelope


def load_scale_band(process, frag="HymnD"):
    """Minimum and maximum over the 7 (mu_F, mu_R) points (see _scale_combos).
    Returns {y: (pt_values, lower, upper)}, or {} if the two mu_F folders are not both there."""
    variation_dirs = [f"{SCALE_DIR}/factor_{factor}" for factor in EXPECTED_SCALE_FACTORS]
    if not all(os.path.isdir(d) for d in variation_dirs):
        return {}

    results_by_combo = {}
    for mu_f, mu_r in _scale_combos():
        results = load_results(process, frag, data_dir=_mu_f_dir(mu_f), mu_r_factor=mu_r)
        if not results:
            return {}
        results_by_combo[(mu_f, mu_r)] = results

    return _envelope_from_results(results_by_combo)


def load_scale_sum_band(frag="HymnD"):
    """Scale-variation band for the diffractive+exclusive sum.
    {} if the two mu_F folders are not both there."""
    variation_dirs = [f"{SCALE_DIR}/factor_{factor}" for factor in EXPECTED_SCALE_FACTORS]
    if not all(os.path.isdir(d) for d in variation_dirs):
        return {}

    totals_by_combo = {}
    for mu_f, mu_r in _scale_combos():
        totals = _summed_results(_mu_f_dir(mu_f), frag, mu_r_factor=mu_r)
        if not totals:
            return {}
        totals_by_combo[(mu_f, mu_r)] = totals

    return _envelope_from_results(totals_by_combo)


def _combine_in_quadrature(central, rep_band, bk_band, scale_band):
    """Total uncertainty around the central curve: sqrt of the sum of squares (now
    only the scale variation). Returns {y: (pt_values, central_values, sigma_total)}."""
    if not central or not scale_band:
        return {}

    combined = {}
    for y in central:
        if y not in scale_band:
            continue
        central_points = dict(central[y])
        # rep_pt, _, rep_stds = rep_band[y]
        # rep_lookup = dict(zip(rep_pt, rep_stds))
        # bk_pt, _, bk_stds = bk_band[y]
        # bk_lookup = dict(zip(bk_pt, bk_stds))
        scale_pt, scale_lower, scale_upper = scale_band[y]
        scale_lookup = {pt: (upper - lower) / 2 for pt, lower, upper in zip(scale_pt, scale_lower, scale_upper)}

        pt_values = sorted(set(central_points) & set(scale_lookup))
        central_values = [central_points[pt] for pt in pt_values]
        sigma_total = [
            # math.sqrt(rep_lookup[pt]**2 + bk_lookup[pt]**2 + scale_lookup[pt]**2)  # replica + BK-IC + scale
            scale_lookup[pt]
            for pt in pt_values
        ]
        combined[y] = (pt_values, central_values, sigma_total)

    return combined


def load_hymnd_band(process, frag="HymnD"):
    """Total uncertainty band (now only the scale variation). {} if there is no data."""
    return _combine_in_quadrature(
        load_results(process, frag),
        load_hymnd_replica_band(process, frag),
        load_bk_band(process, frag),
        load_scale_band(process, frag),
    )


def load_hymnd_sum_band(frag="HymnD"):
    """"HymnD" total band (see load_hymnd_band) for the diffractive+
    exclusive SUM.
    """
    return _combine_in_quadrature(
        _summed_results(CENTRAL_DIR, frag),
        load_hymnd_replica_sum_band(frag),
        load_bk_sum_band(frag),
        load_scale_sum_band(frag),
    )


def main():
    all_results = {}
    for process in PROCESSES:
        for frag in FRAG_TYPES:
            results = load_results(process, frag)
            if results:
                all_results[(process, frag)] = results

    hymnd_bands = {process: load_hymnd_band(process) for process in PROCESSES}

    if not all_results and not any(hymnd_bands.values()):
        sys.exit(f"No files found for NUCLEUS={NUCLEUS}, CHANNEL={CHANNEL} "
                  "in ../output/ -- run ../run_scripts/run_nucleus.sh first.")

    def y_in_range(y):
        return MIN_Y_TO_PLOT <= y <= MAX_Y_TO_PLOT and y == round(y)

    y_values = set()
    for results in all_results.values():
        for y in results:
            if y_in_range(y):
                y_values.add(y)
    for band in hymnd_bands.values():
        for y in band:
            if y_in_range(y):
                y_values.add(y)
    y_values = sorted(y_values)

    palette = ["#053061", "#2166ac", "#67a9cf", "#ef8a62", "#b2182b"]
    y_colors = dict(zip([-2.0, -1.0, 0.0, 1.0, 2.0], palette))
    colors = {y: y_colors.get(y, palette[i % len(palette)]) for i, y in enumerate(y_values)}

    plt.figure(figsize=(8, 7))

    for process, frag in LINESTYLES:
        linestyle = LINESTYLES[(process, frag)]

        hymnd_band = hymnd_bands.get(process) if frag == "HymnD" else None
        if hymnd_band:
            for y in sorted(hymnd_band):
                if not y_in_range(y):
                    continue
                pt_values, central_values, sigma_total = hymnd_band[y]
                lower = [v - s for v, s in zip(central_values, sigma_total)]
                upper = [v + s for v, s in zip(central_values, sigma_total)]
                plt.fill_between(pt_values, lower, upper, color=colors[y], alpha=0.25, linewidth=0)
                plt.plot(pt_values, central_values, color=colors[y], linestyle=linestyle, lw=LINEWIDTH)
            continue

        results = all_results.get((process, frag))
        if not results:
            continue
        for y in sorted(results):
            if not y_in_range(y):
                continue
            points = sorted(results[y])
            pt_values = [pair[0] for pair in points]
            cross_section_values = [pair[1] for pair in points]
            plt.plot(pt_values, cross_section_values, color=colors[y], linestyle=linestyle, lw=LINEWIDTH)

    plt.yscale("log")
    plt.ylim(3e-9, 1.5 * plt.gca().dataLim.y1)   # just above the highest curve
    plt.xlim(0, 12)
    plt.xlabel(r"$p_{D^0\perp}$ [GeV]", labelpad=14)
    plt.ylabel(r"$d\sigma/dy\,dp_{D^0\perp}$ [mb/GeV]", labelpad=16)
    plt.text(0.95, 0.92, f"{NUCLEUS}-{NUCLEUS} 5.36 TeV\n{CHANNEL}, HymnD", transform=plt.gca().transAxes,
             ha="right", va="top", fontsize=24, linespacing=1.8)
    plt.gca().yaxis.set_minor_locator(
        ticker.LogLocator(base=10.0, subs=[2, 3, 4, 5, 6, 7, 8, 9], numticks=100))
    plt.gca().yaxis.set_minor_formatter(ticker.NullFormatter())
    plt.gca().yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=100))
    plt.gca().yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda val, pos: f"$10^{{{round(math.log10(val))}}}$"
                              if round(math.log10(val)) % 2 != 0 else ""))

    y_handles = [Line2D([0], [0], color=colors[y], linestyle="-", linewidth=3, label=f"$y={y:g}$")
                 for y in y_values]
    y_legend = plt.legend(handles=y_handles, loc="lower left", bbox_to_anchor=(0.02, 0.0), fontsize=22,
                          frameon=False)
    plt.gca().add_artist(y_legend)

    style_handles = []
    process_labels = {"exclusive": "Exclusive", "diffractive": r"Diffractive$_{\mbox{\fontsize{13.5}{13.5}\selectfont TMD}}$"}
    for process in ["exclusive", "diffractive"]:
        for frag in ["HymnD", "BCFY", "KniehlKramer"]:
            if (process, frag) not in LINESTYLES:
                continue
            linestyle = LINESTYLES[(process, frag)]
            label = process_labels[process] if frag == "HymnD" else f"{process_labels[process]}, {frag}"
            if frag == "HymnD" and hymnd_bands.get(process):
                style_handles.append(Line2D([0], [0], color="0.3", linestyle=linestyle, lw=LINEWIDTH, label=label))
            elif (process, frag) in all_results:
                style_handles.append(Line2D([0], [0], color="0.3", linestyle=linestyle, lw=LINEWIDTH, label=label))
    # under the text in the upper right
    plt.legend(handles=style_handles, loc="upper right", bbox_to_anchor=(0.98, 0.76), fontsize=22,
               frameon=False)

    plt.tight_layout()
    outname = f"../plots/D0_{CHANNEL}_{NUCLEUS}{NUCLEUS}.pdf"
    plt.savefig(outname, bbox_inches="tight")
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
