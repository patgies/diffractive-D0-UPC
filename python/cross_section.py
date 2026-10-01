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
sys.path.insert(0, os.path.dirname(__file__))
from alphas_running import alphas_run

# make the plot look nicer (same style as xpom_spectrum.py / fixed_qp_spectrum.py)
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

# This script reads the data files written by ../local_workflows/run_many_nucleus.sh
# (produced by the D0 program), integrates them over the target-nucleus
# impact parameter b using Simpson's rule, and makes one plot with every
# combination of process (exclusive/diffractive) and fragmentation function
# (BCFY/KniehlKramer) on it. Each rapidity gets its own color, and each
# (process, fragmentation) combination gets its own linestyle.

alphae = 1/137
mc     = 1.5          # charm mass in GeV
e_c    = 2/3
Nc     = 3

# conversion factor from GeV^-2 to mb
FMGEV = 5.068
GEVSQR_TO_NB = 1.0e7 / (FMGEV * FMGEV)
GEVSQR_TO_MB = GEVSQR_TO_NB * 1e-6

NUCLEUS = os.environ.get("NUCLEUS", "Pb")
CHANNEL = os.environ.get("CHANNEL", "An0n").translate(str.maketrans('', '', '() '))

PROCESSES = ["diffractive", "exclusive"]
FRAG_TYPES = ["HymnD"]
MIN_Y_TO_PLOT = -1.0  # keeps the plot from getting too crowded
MAX_Y_TO_PLOT = 3.0   # keeps the plot from getting too crowded

# mu_F (fragmentation/factorization scale, Q = factor*mt0 -- see SCALE_FACTOR
# in src/main.cpp/main_xpom.cpp) and mu_R (renormalization scale used inside
# alpha_s, see prefactor() below) each range over these three factors times
# mu0 = mT = sqrt(pD0^2+mc^2), following e.g. Cacciari, Innocenti & Stasto,
# arXiv:2506.09893. mu_F=1.0 tags to the plain (no-subdir) run; the two
# non-central mu_F values need their own pre-computed run_many_nucleus.sh
# output (run_HymnD_scale_variation.sh), since Q feeds into the HymnD FF
# grid read by the C++ side. mu_R has no such directory: alpha_s(mu_R) is
# applied here in Python as a plain multiplicative factor in prefactor(), on
# top of whichever mu_F directory's raw dsigma is being read -- so varying
# mu_R alone never needs a new VEGAS run.
SCALE_FACTOR_TAGS = {0.5: "0.5", 1.0: None, 2.0: "2.0"}
SCALE_FACTORS = list(SCALE_FACTOR_TAGS)                  # [0.5, 1.0, 2.0]
EXPECTED_SCALE_FACTORS = [tag for tag in SCALE_FACTOR_TAGS.values() if tag]  # ["0.5", "2.0"]
# must match run_HymnD_scale_variation.sh's SCALE_FACTORS


def _mu_f_dir(factor):
    """Raw-output directory for a given mu_F factor (see SCALE_FACTOR_TAGS)."""
    tag = SCALE_FACTOR_TAGS[factor]
    return ".." if tag is None else f"../files/HymnD_scale/factor_{tag}"

# which linestyle to use for each (process, fragmentation) combination
LINESTYLES = {
    ("diffractive", "BCFY"):          "-",
    ("diffractive", "KniehlKramer"):  "--",
    ("diffractive", "HymnD"):        "-",
    ("exclusive",   "BCFY"):          ":",
    ("exclusive",   "KniehlKramer"):  "-.",
    ("exclusive",   "HymnD"):        ":",
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


def prefactor(process, pt, mu_r_factor=1.0):
    """Physical prefactor, matching plot_pt_spectrum.py's convention.
    No sigma0 here: that's the GBW proton-dipole normalization, and for a
    nucleus target the "area" is already accounted for by integrate_over_b's
    Glauber b-integral -- applying sigma0 on top of that would double-count it.

    mu_r_factor rescales the renormalization scale mu_R = mu_r_factor*mT fed
    into alpha_s -- independent of the mu_F used for the fragmentation
    function (that one lives in which raw-output directory load_results()
    reads from, see _mu_f_dir()). Only the diffractive process depends on
    alpha_s, so mu_r_factor is a no-op for exclusive.
    """
    if process == "exclusive":
        return alphae * Nc * e_c**2 / (2 * math.pi**2)
    elif process == "diffractive":
        alphas = alphas_run(mu_r_factor * math.sqrt(pt**2 + mc**2))
        return alphas * alphae * e_c**2 * (Nc**2 - 1) / (8 * math.pi**4)
    raise ValueError(f"unknown process {process}")


def load_results(process, frag, base_dir="..", mu_r_factor=1.0):
    """Read all files for one (process, frag) combination and return {y: [(pt, cross_section), ...]}."""
    pattern = f"{base_dir}/files/D0_{process}_{frag}_{CHANNEL}_{NUCLEUS}_y*.dat"
    results = {}
    for filename in sorted(glob.glob(pattern)):
        y = read_rapidity(filename)
        b_list, pt_list, dsigma_list = read_data_file(filename)
        pt_groups = group_by_pt(b_list, pt_list, dsigma_list)

        # A pT group with fewer b-samples than the file's own full count is
        # a sweep still in progress (e.g. this file was read mid-write by a
        # background run_many_nucleus.sh job) -- integrating over a partial
        # b-grid gives a silently wrong number, not just a noisy one, so
        # skip those points rather than including them. The MIN_B_SAMPLES
        # floor additionally catches the case where EVERY pT group in the
        # file is equally (and uniformly) incomplete -- e.g. a sweep that's
        # only just started, where every point so far has just 1 b-sample,
        # so "matches this file's own max" alone wouldn't flag anything.
        MIN_B_SAMPLES = 10   # real Glauber b-grids here have 17 (Pb/mve) or 18 (bk_posterior) samples
        full_b_count = max((len(pairs) for pairs in pt_groups.values()), default=0)
        pt_groups = {
            pt: pairs for pt, pairs in pt_groups.items()
            if len(pairs) == full_b_count and len(pairs) >= MIN_B_SAMPLES
        }

        if not pt_groups:
            continue   # nothing survived the completeness filter for this y -- don't add an empty entry

        points = []
        for pt, pairs in sorted(pt_groups.items()):
            b_integral = integrate_over_b(pairs)
            # 2*pi*b_integral is the Glauber transverse-plane (b) integral.
            # 2*pi*pt is the Jacobian from d^2pD0 to dpD0.
            cross_section = (2*math.pi) * b_integral * prefactor(process, pt, mu_r_factor) \
                            * (2*math.pi) * pt * GEVSQR_TO_MB
            points.append((pt, cross_section))
        results[y] = points
    return results


def _band_from_members(process, frag, member_dir_pattern):
    """Shared aggregation for load_hymnd_replica_band/load_bk_band: combine one
    run_many_nucleus.sh output set per member directory into a mean +/-
    standard-deviation band per rapidity (both uncertainty sources here are
    sampled sets -- HymnD replicas, BK posterior samples -- not Hessian
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


def load_hymnd_replica_band(process, frag="HymnD"):
    """HymnD-fragmentation-replica band (see run_HymnD_members_roihu.sbatch):
    dipole amplitude held fixed at the central data/Pb/mve/ set, fragmentation
    function varied across the 101 replica members. {} until
    run_HymnD_members_roihu.sbatch's output has been copied back into
    ../files/HymnD/member_*.
    """
    return _band_from_members(process, frag, "../files/HymnD/member_*")


def load_bk_band(process, frag="HymnD"):
    """BK-initial-condition posterior band (see
    run_bk_posterior_members_roihu.sbatch): fragmentation function held fixed
    at the central HymnD member, dipole amplitude varied across the 100
    Bayesian BK-IC posterior samples in bk/. {} until
    run_bk_posterior_members_roihu.sbatch's output has been copied back into
    ../files/bk_posterior/member_*.

    Currently unused in load_hymnd_band/load_hymnd_sum_band -- see the
    commented-out lines there.
    """
    return _band_from_members(process, frag, "../files/bk_posterior/member_*")


def _summed_results(base_dir, frag="HymnD", mu_r_factor=1.0):
    """diffractive+exclusive summed per (y, pt) for one directory. Returns
    {y: [(pt, total), ...]}, or {} if neither process has data there.
    """
    diff = load_results("diffractive", frag, base_dir=base_dir, mu_r_factor=mu_r_factor)
    excl = load_results("exclusive", frag, base_dir=base_dir, mu_r_factor=mu_r_factor)
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
    """Like _band_from_members, but sums diffractive+exclusive PER MEMBER
    (via _summed_results) before aggregating -- statistically correct,
    since both processes in a given member share the same FF replica /
    BK-IC dipole sample (they're not independent draws), unlike summing two
    independently-aggregated per-process bands would be.

    Returns {y: (pt_values, means, stds)}, or {} if no member directories
    are found yet.
    """
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
    """HymnD-replica band for the diffractive+exclusive SUM (see
    _summed_band_from_members). {} until run_HymnD_members_roihu.sbatch's
    output is in ../files/HymnD/member_*.
    """
    return _summed_band_from_members("../files/HymnD/member_*", frag)


def load_bk_sum_band(frag="HymnD"):
    """BK-IC posterior band for the diffractive+exclusive SUM (see
    _summed_band_from_members). {} until run_bk_posterior_members_roihu.sbatch's
    output is in ../files/bk_posterior/member_*.

    Currently unused in load_hymnd_sum_band -- see the commented-out lines
    there.
    """
    return _summed_band_from_members("../files/bk_posterior/member_*", frag)


def _scale_combos():
    """The restricted 7-point (mu_F, mu_R) grid: both range over
    SCALE_FACTORS, dropping the two anti-correlated corners mu_F/mu_R in
    {4, 1/4} -- the conventional 7-point scheme (e.g. arXiv:2506.09893),
    since those two extreme combinations aren't believed to give a
    meaningful estimate of missing higher orders.
    """
    return [(mu_f, mu_r) for mu_f in SCALE_FACTORS for mu_r in SCALE_FACTORS
            if 0.5 <= mu_f / mu_r <= 2.0]


def _envelope_from_results(results_by_combo):
    """Shared min/max-over-combos envelope builder for load_scale_band/
    load_scale_sum_band. results_by_combo: {(mu_f, mu_r): {y: [(pt, val), ...]}}.
    Returns {y: (pt_values, lower, upper)}.
    """
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
    """Combined factorization/renormalization-scale-variation envelope: the
    restricted 7-point (mu_F, mu_R) grid (see _scale_combos), mu_F, mu_R in
    {0.5, 1, 2} * mu0 (mu0 = mT = sqrt(pD0^2+mc^2)), following e.g.
    Cacciari, Innocenti & Stasto, arXiv:2506.09893.

    mu_F only affects the HymnD fragmentation function's evolution scale Q,
    which is baked into the raw dsigma files at generation time (see
    SCALE_FACTOR in src/main.cpp/main_xpom.cpp, run via
    run_HymnD_scale_variation.sh) -- so its two non-central values need
    their own pre-computed run_many_nucleus.sh output directories
    (_mu_f_dir). mu_R only enters through alpha_s(mu_R) in prefactor(),
    applied here in Python on top of whichever mu_F directory's raw dsigma
    is being read -- so varying mu_R needs no extra VEGAS run, just a
    different mu_r_factor passed to load_results() on the SAME 3 files.

    This is a scale-convention envelope (min/max over the 7 points), not a
    statistical replica sample, so it's built differently from
    load_hymnd_replica_band -- min and max, not mean +/- std.

    Returns {y: (pt_values, lower, upper)}, or {} until BOTH non-central
    mu_F directories exist in ../files/HymnD_scale/ (see
    EXPECTED_SCALE_FACTORS -- checking each by name rather than just
    globbing "factor_*" and using whatever's there matters here:
    run_HymnD_scale_variation.sh runs the two factors sequentially, so
    factor_2.0's directory doesn't exist on disk AT ALL until factor_0.5 has
    fully finished -- glob-and-use-whatever-exists would silently treat a
    still-running factor_0.5 alone as if it were the complete envelope).
    """
    variation_dirs = [f"../files/HymnD_scale/factor_{factor}" for factor in EXPECTED_SCALE_FACTORS]
    if not all(os.path.isdir(d) for d in variation_dirs):
        return {}

    results_by_combo = {}
    for mu_f, mu_r in _scale_combos():
        results = load_results(process, frag, base_dir=_mu_f_dir(mu_f), mu_r_factor=mu_r)
        if not results:
            return {}
        results_by_combo[(mu_f, mu_r)] = results

    return _envelope_from_results(results_by_combo)


def load_scale_sum_band(frag="HymnD"):
    """Scale-variation envelope (see load_scale_band) for the diffractive+
    exclusive SUM. {} until BOTH non-central mu_F directories' output exists
    (see load_scale_band's docstring for why this checks by name rather
    than globbing).
    """
    variation_dirs = [f"../files/HymnD_scale/factor_{factor}" for factor in EXPECTED_SCALE_FACTORS]
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
    """Shared combination for load_hymnd_band/load_hymnd_sum_band: adds
    uncertainty sources in quadrature around the plain central curve.
    Currently only the scale-variation envelope half-width is active
    (sigma_total = sigma_scale) -- the HymnD-replica std and BK-IC
    uncertainty are both commented out of the sum below, but rep_band/
    bk_band are still threaded through so re-enabling either is a
    one-line change. Requires the (currently) active source present
    (returns {} otherwise), since a partial combination would understate
    the true uncertainty.

    Returns {y: (pt_values, central_values, sigma_total)}.
    """
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
    """"HymnD" combined theory-uncertainty band: currently just the
    fragmentation-scale-variation uncertainty (see _combine_in_quadrature
    -- replica and BK-IC are both commented out there). {} until
    load_scale_band has data.
    """
    return _combine_in_quadrature(
        load_results(process, frag),
        load_hymnd_replica_band(process, frag),
        load_bk_band(process, frag),
        load_scale_band(process, frag),
    )


def load_hymnd_sum_band(frag="HymnD"):
    """"HymnD" combined band (see load_hymnd_band) for the diffractive+
    exclusive SUM.
    """
    return _combine_in_quadrature(
        _summed_results("..", frag),
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

    # If run_HymnD_members_roihu.sbatch's / run_HymnD_scale_variation.sh's
    # output are both present, plot the combined "HymnD" theory-uncertainty
    # band (HymnD-replica and scale-variation uncertainties added in
    # quadrature -- see load_hymnd_band) instead of the single-member
    # HymnD line.
    hymnd_bands = {process: load_hymnd_band(process) for process in PROCESSES}

    if not all_results and not any(hymnd_bands.values()):
        sys.exit(f"No files found for NUCLEUS={NUCLEUS}, CHANNEL={CHANNEL} "
                  "in ../files/ -- run ../local_workflows/run_many_nucleus.sh first.")

    # collect every rapidity value that shows up in any of the results,
    # keeping only integer y in [MIN_Y_TO_PLOT, MAX_Y_TO_PLOT] so the plot
    # doesn't get too crowded
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

        hymnd_band = hymnd_bands.get(process) if frag == "HymnD" else None
        if hymnd_band:
            for y in sorted(hymnd_band):
                if not y_in_range(y):
                    continue
                pt_values, central_values, sigma_total = hymnd_band[y]
                lower = [v - s for v, s in zip(central_values, sigma_total)]
                upper = [v + s for v, s in zip(central_values, sigma_total)]
                plt.fill_between(pt_values, lower, upper, color=colors[y], alpha=0.25, linewidth=0)
                plt.plot(pt_values, central_values, color=colors[y], linestyle=linestyle)
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
            plt.plot(pt_values, cross_section_values, color=colors[y], linestyle=linestyle)

    plt.yscale("log")
    # Fixed range rather than autoscale: at very low pT, mt0~m_c and
    # SCALE_FACTOR=0.5 pushes Q below the charm-mass threshold, where the
    # fragmentation function genuinely goes to ~0 (real physics, see
    # notes/) -- that blows up the HymnD band's lower edge toward the
    # 1e-30 floor at a few low-pT points, which would otherwise stretch
    # the whole axis and make everything else unreadable.
    plt.ylim(1e-10, 1e2)
    plt.xlabel(r"$p_{D^0\perp}$ [GeV]", labelpad=15)
    plt.ylabel(r"$d\sigma/dy\,dp_{D^0\perp}$ [mb/GeV]", labelpad=15)
    plt.title(f"$D^0$ photoproduction, {NUCLEUS}+{NUCLEUS} UPC ({CHANNEL}, HymnD)", pad=15)
    # Force log-scale minor ticks (2,3,...,9 within each decade) to actually
    # show. LogLocator has an internal numticks budget that silently
    # returns NO minor ticks at all once the axis spans too many decades
    # (ours spans ~10) -- numticks=100 raises that budget so it doesn't
    # give up (default budget is far too small for this range).
    plt.gca().yaxis.set_minor_locator(
        ticker.LogLocator(base=10.0, subs=[2, 3, 4, 5, 6, 7, 8, 9], numticks=100))
    plt.gca().yaxis.set_minor_formatter(ticker.NullFormatter())
    # Same numticks issue for the major locator: with ~10 decades, the
    # default budget skips tick MARKS (not just labels) at every other
    # decade. Force a tick at every decade, full major-tick length, but
    # only label every other one (matches the space we actually have).
    plt.gca().yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=100))
    plt.gca().yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda val, pos: f"$10^{{{round(math.log10(val))}}}$"
                              if round(math.log10(val)) % 2 == 0 else ""))

    # Discrete per-y legend (a plain colorbar doesn't read well once y is
    # restricted to a handful of integers) -- a separate legend box from
    # the process/frag one below, added via add_artist so the second
    # plt.legend() call doesn't replace it.
    y_handles = [Line2D([0], [0], color=colors[y], linestyle="-", linewidth=2, label=f"$y={y:g}$")
                 for y in y_values]
    y_legend = plt.legend(handles=y_handles, loc="upper right", fontsize=15)
    plt.gca().add_artist(y_legend)

    style_handles = []
    for process in ["exclusive", "diffractive"]:
        for frag in ["HymnD", "BCFY", "KniehlKramer"]:
            if (process, frag) not in LINESTYLES:
                continue
            linestyle = LINESTYLES[(process, frag)]
            label = process if frag == "HymnD" else f"{process}, {frag}"
            if frag == "HymnD" and hymnd_bands.get(process):
                style_handles.append(Line2D([0], [0], color="black", linestyle=linestyle, label=label))
            elif (process, frag) in all_results:
                style_handles.append(Line2D([0], [0], color="black", linestyle=linestyle, label=label))
    plt.legend(handles=style_handles, loc="lower left", fontsize=15)

    plt.tight_layout()
    outname = f"../plots/cross_section_{CHANNEL}_{NUCLEUS}.pdf"
    plt.savefig(outname, bbox_inches="tight")
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
