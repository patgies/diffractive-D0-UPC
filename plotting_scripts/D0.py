import glob
import math
import os
import statistics
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
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
    "xtick.color": "gray",
    "ytick.color": "gray",
    "xtick.labelcolor": "black",
    "ytick.labelcolor": "black",
    "xtick.major.pad": 8,
    "ytick.major.pad": 4,
})

# Reads the D0 output of run_nucleus.sh, integrates it over the impact parameter b_d
# and plots exclusive and diffractive spectra. The other plotting scripts take the plot style,
# the constants, the prefactors and the readers from here.

alphae = 1/137
mc     = 1.5
e_c    = 2/3
Nc     = 3
SIGMA0_MB = 16.36   # sigma0/2 of the proton dipole fit (MVe), in mb: the proton transverse area
A_NUCLEUS = {"Pb": 208}

FMGEV = 5.068
GEVSQR_TO_NB = 1.0e7 / (FMGEV * FMGEV)
GEVSQR_TO_MB = GEVSQR_TO_NB * 1e-6

NUCLEUS = os.environ.get("NUCLEUS", "Pb")
CHANNEL = os.environ.get("CHANNEL", "An0n").translate(str.maketrans('', '', '() '))

# FLUX_MODEL: EFF (default) or another flux, whose files are in ../output/<FLUX_MODEL>/ with a _<FLUX_MODEL> tag
FLUX_MODEL = os.environ.get("FLUX_MODEL", "EFF")
FLUX_TAG = "" if FLUX_MODEL == "EFF" else f"_{FLUX_MODEL}"
OUTPUT_BASE = "../output" if FLUX_MODEL == "EFF" else f"../output/{FLUX_MODEL}"

CENTRAL_DIR = os.environ.get("CENTRAL_DIR", f"{OUTPUT_BASE}/{CHANNEL}/central_values")
SCALE_DIR = os.environ.get("SCALE_DIR", f"{OUTPUT_BASE}/{CHANNEL}/scale_variation/HymnD")

PROCESSES = ["diffractive", "exclusive"]
FRAG_TYPES = ["HymnD"]
MIN_Y_TO_PLOT = -1.0
MAX_Y_TO_PLOT = 2.0
MIN_B_D_SAMPLES = 10   # a pT point needs at least this many impact parameters to be integrated

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


# Dilute edge of the nucleus. The dipole files stop at the impact parameter b_max where
# A T_A(b) sigma0/2 ~ 1; beyond it the nuclear amplitude is A T_A(b) sigma0/2 times the proton one
# (arXiv:1309.6963), so that for a coherent process sigma_A^edge = nuclear_edge_factor * sigma_p.
# sigma_p: proton target with the same flux, channel, fragmentation function and scale, in
# ../output/proton_edge/<FLUX_MODEL>/<CHANNEL>/factor_<scale>/<FRAG>/ (run_proton.sh, TARGET=AA).
NUCLEAR_EDGE = os.environ.get("NUCLEAR_EDGE", "1") == "1"
PROTON_EDGE_DIR = os.environ.get("PROTON_EDGE_DIR", "../output/proton_edge")
WOODS_SAXON = {"Pb": (6.49, 0.542)}   # R and a in fm that reproduce the initial conditions of the dipole files
HBARC = 0.1973   # GeV fm
_edge_cache = {}


def nuclear_edge_factor(b_max, sigma0_mb=SIGMA0_MB, nucleus=None):
    """(sigma0/2) * integral of [A T_A(b)]^2 d^2b over b > b_max (b_max in GeV^-1)."""
    nucleus = NUCLEUS if nucleus is None else nucleus
    key = (nucleus, round(b_max, 6))
    if key not in _edge_cache:
        radius, diffuseness = WOODS_SAXON[nucleus]
        z = np.linspace(0.0, 30.0, 3001)
        b = np.linspace(0.0, 25.0, 2501)
        thickness = np.array([2 * np.trapz(1 / (1 + np.exp((np.sqrt(bi**2 + z**2) - radius) / diffuseness)), z)
                              for bi in b])
        thickness *= A_NUCLEUS[nucleus] / np.trapz(2 * math.pi * b * thickness, b)   # A T_A in fm^-2
        edge = b >= b_max * HBARC
        _edge_cache[key] = np.trapz(2 * math.pi * b[edge] * thickness[edge]**2, b[edge])   # fm^-2
    return sigma0_mb / 10.0 * _edge_cache[key]   # 1 fm^2 = 10 mb


def proton_prefactor(process, pt, mu_r_factor=1.0):
    """Prefactor for a proton target, with sigma0/2 in mb (so no GEVSQR_TO_MB)."""
    if process == "exclusive":
        return alphae * Nc * e_c**2 * SIGMA0_MB / (2 * math.pi**2)
    alphas = alphas_run(mu_r_factor * math.sqrt(pt**2 + mc**2))
    return alphas * alphae * e_c**2 * (Nc**2 - 1) * SIGMA0_MB / (8 * math.pi**4)


_warned_no_edge = set()


def load_proton_edge(process, frag, data_dir, mu_r_factor=1.0):
    """Proton-target cross section {y: {pt: dsigma/dy dpt in mb/GeV}} that matches the lead results of
    data_dir: the scale factor is read from a "factor_<f>" folder name (1.0 otherwise; the member
    folders of the HymnD and BK bands use the central proton result). {} if the runs are missing."""
    scale = next((part[len("factor_"):] for part in data_dir.split("/") if part.startswith("factor_")), "1.0")
    folder = f"{PROTON_EDGE_DIR}/{FLUX_MODEL}/{CHANNEL}/factor_{scale}/{frag}"
    files = sorted(glob.glob(f"{folder}/D0_proton_baseline_{process}_{frag}_{CHANNEL}_y*.dat"))
    if not files:
        if folder not in _warned_no_edge:
            _warned_no_edge.add(folder)
            print(f"Warning: no proton results in {folder} -- the nuclear edge is not added.", file=sys.stderr)
        return {}
    results = {}
    for filename in files:
        points = {}
        with open(filename) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    pt, raw = (float(v) for v in line.split())
                    points[pt] = raw * proton_prefactor(process, pt, mu_r_factor) * (2 * math.pi) * pt
        results[read_rapidity(filename)] = points
    return results


def load_results(process, frag, data_dir=None, mu_r_factor=1.0):
    """Read all files for one (process, frag) in data_dir or data_dir/<frag>/
    and return {y: [(pt, cross_section), ...]}."""
    data_dir = CENTRAL_DIR if data_dir is None else data_dir
    name = f"D0_{process}_{frag}_{CHANNEL}_{NUCLEUS}{FLUX_TAG}_y*.dat"
    filenames = glob.glob(f"{data_dir}/{name}") or glob.glob(f"{data_dir}/{frag}/{name}")
    results = {}
    add_edge = NUCLEAR_EDGE and NUCLEUS in A_NUCLEUS and bool(filenames)
    proton_edge = load_proton_edge(process, frag, data_dir, mu_r_factor) if add_edge else {}
    for filename in sorted(filenames):
        y = read_rapidity(filename)
        b_d_list, pt_list, dsigma_list = read_data_file(filename)
        pt_groups = group_by_pt(b_d_list, pt_list, dsigma_list)

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
            if pt in proton_edge.get(y, {}):
                cross_section += nuclear_edge_factor(max(b_d for b_d, _ in pairs)) * proton_edge[y][pt]
            elif proton_edge and filename not in _warned_no_edge:   # the proton run lacks this y or pt
                _warned_no_edge.add(filename)
                print(f"Warning: no proton result for y={y:g}, pt={pt:g} ({filename}) -- "
                      "the nuclear edge is not added there.", file=sys.stderr)
            points.append((pt, cross_section))
        results[y] = points
    return results


def _member_band(member_dir_pattern, load):
    """Mean and standard deviation over the member folders of load(member_dir), which gives
    {y: [(pt, value), ...]}. Returns {y: (pt_values, means, stds)}, or {} if none are found."""
    per_member = [load(member_dir) for member_dir in sorted(glob.glob(member_dir_pattern))]
    per_member = [results for results in per_member if results]
    if not per_member:
        return {}

    band = {}
    for y in per_member[0]:
        pt_values = sorted(pt for pt, _ in per_member[0][y])
        values_by_pt = [[dict(results[y])[pt] for results in per_member if y in results] for pt in pt_values]
        means = [statistics.mean(values) for values in values_by_pt]
        stds = [statistics.pstdev(values) for values in values_by_pt]
        band[y] = (pt_values, means, stds)
    return band


def _band_from_members(process, frag, member_dir_pattern):
    """Band of one process from the member folders (see _member_band)."""
    return _member_band(member_dir_pattern, lambda member_dir: load_results(process, frag, data_dir=member_dir))


def load_hymnd_replica_band(process, frag="HymnD"):
    """Band from the HymnD replicas (one fragmentation function per member). {} if there is no data."""
    return _band_from_members(process, frag, f"{OUTPUT_BASE}/{CHANNEL}/HymnD_band/member_*")


def load_bk_band(process, frag="HymnD"):
    """Band from the BK posterior samples (one dipole per sample). {} if there is no data."""
    return _band_from_members(process, frag, f"{OUTPUT_BASE}/{CHANNEL}/bk_band/member_*")


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
    """Like _band_from_members, but for the diffractive+exclusive sum of each member."""
    return _member_band(member_dir_pattern, lambda member_dir: _summed_results(member_dir, frag))


def load_hymnd_replica_sum_band(frag="HymnD"):
    """Band from the HymnD replicas for the diffractive+exclusive sum. {} if there is no data."""
    return _summed_band_from_members(f"{OUTPUT_BASE}/{CHANNEL}/HymnD_band/member_*", frag)


def load_bk_sum_band(frag="HymnD"):
    """Band from the BK posterior samples for the diffractive+exclusive sum. {} if there is no data."""
    return _summed_band_from_members(f"{OUTPUT_BASE}/{CHANNEL}/bk_band/member_*", frag)


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


def _scale_envelope(load):
    """Minimum and maximum over the 7 (mu_F, mu_R) points (see _scale_combos) of
    load(data_dir, mu_r_factor) = {y: [(pt, value), ...]}. Returns {y: (pt_values, lower, upper)},
    or {} if the two mu_F folders are not both there."""
    variation_dirs = [f"{SCALE_DIR}/factor_{factor}" for factor in EXPECTED_SCALE_FACTORS]
    if not all(os.path.isdir(d) for d in variation_dirs):
        return {}

    results_by_combo = {}
    for mu_f, mu_r in _scale_combos():
        results = load(_mu_f_dir(mu_f), mu_r)
        if not results:
            return {}
        results_by_combo[(mu_f, mu_r)] = results

    return _envelope_from_results(results_by_combo)


def load_scale_band(process, frag="HymnD"):
    """Scale-variation band of one process (see _scale_envelope)."""
    return _scale_envelope(lambda data_dir, mu_r: load_results(process, frag, data_dir=data_dir, mu_r_factor=mu_r))


def load_scale_sum_band(frag="HymnD"):
    """Scale-variation band for the diffractive+exclusive sum (see _scale_envelope)."""
    return _scale_envelope(lambda data_dir, mu_r: _summed_results(data_dir, frag, mu_r_factor=mu_r))


def _relative_std(member_band, y):
    """{pt: std/mean} of a member band (see _band_from_members) at rapidity y; {} if it has no data."""
    if not member_band or y not in member_band:
        return {}
    pt_values, means, stds = member_band[y]
    return {pt: std / mean for pt, mean, std in zip(pt_values, means, stds) if mean > 0}


def _combine_in_quadrature(central, rep_band, bk_band, scale_band):
    """Total uncertainty around the central curve: scale variation, HymnD replicas and BK posterior
    samples added in quadrature. The replica and BK bands enter as relative spreads (std/mean of
    the members), since the BK samples come from another fit than the central dipole; a band with
    no data (other channels) is left out. Returns {y: (pt_values, central_values, sigma_total)}."""
    if not central or not scale_band:
        return {}

    combined = {}
    for y in central:
        if y not in scale_band:
            continue
        central_points = dict(central[y])
        rep_lookup = _relative_std(rep_band, y)
        bk_lookup = _relative_std(bk_band, y)
        scale_pt, scale_lower, scale_upper = scale_band[y]
        scale_lookup = {pt: (upper - lower) / 2 for pt, lower, upper in zip(scale_pt, scale_lower, scale_upper)}

        pt_values = sorted(set(central_points) & set(scale_lookup))
        central_values = [central_points[pt] for pt in pt_values]
        sigma_total = [
            math.sqrt(scale_lookup[pt]**2
                      + (rep_lookup.get(pt, 0.0) * central_points[pt])**2
                      + (bk_lookup.get(pt, 0.0) * central_points[pt])**2)
            for pt in pt_values
        ]
        combined[y] = (pt_values, central_values, sigma_total)

    return combined


def band_label():
    """Legend label of the uncertainty band: the total one where the HymnD-replica and BK runs exist
    for this channel (see _combine_in_quadrature), the scale variation alone otherwise."""
    has_members = (glob.glob(f"{OUTPUT_BASE}/{CHANNEL}/HymnD_band/member_*")
                   or glob.glob(f"{OUTPUT_BASE}/{CHANNEL}/bk_band/member_*"))
    if has_members:
        return r"Scale $\oplus$ FF $\oplus$ dipole fit"
    return r"$\mu_F, \mu_R \in [0.5, 2]\, m_t$"


def load_hymnd_band(process, frag="HymnD"):
    """Total uncertainty band (scale + HymnD replicas + BK samples, see _combine_in_quadrature).
    {} if there is no data."""
    return _combine_in_quadrature(
        load_results(process, frag),
        load_hymnd_replica_band(process, frag),
        load_bk_band(process, frag),
        load_scale_band(process, frag),
    )


def load_hymnd_sum_band(frag="HymnD"):
    """Total uncertainty band (see load_hymnd_band) for the diffractive+exclusive sum."""
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
    plt.ylim(1e-11, 1.5 * plt.gca().dataLim.y1)   # just above the highest curve; room for the legends below
    plt.xlim(0, 12)
    plt.xlabel(r"$p_{D\perp}$ [GeV]", labelpad=14)
    plt.ylabel(r"$d\sigma/dy\,dp_{D\perp}$ [mb/GeV]", labelpad=11)
    plt.text(0.95, 0.945, f"{NUCLEUS}-{NUCLEUS} 5.36 TeV\n{CHANNEL}, HymnD", transform=plt.gca().transAxes,
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
    if any(hymnd_bands.values()):
        style_handles.append(Patch(facecolor="0.3", alpha=0.25, edgecolor="none", label=band_label()))
    # lower right, below the exclusive curves
    plt.legend(handles=style_handles, loc="lower right", bbox_to_anchor=(1.01, 0.0), fontsize=22,
               frameon=False)

    plt.tight_layout()
    outname = f"../plots/D0_{CHANNEL}_{NUCLEUS}{NUCLEUS}.pdf"
    plt.savefig(outname, bbox_inches="tight")
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
