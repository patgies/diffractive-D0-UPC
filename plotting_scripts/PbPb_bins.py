import os
import sys


from D0 import load_results, _summed_results, _scale_combos, SCALE_FACTOR_TAGS, NUCLEUS, CHANNEL, CENTRAL_DIR
from pPb_bins import bin_averages, draw_panels, FRAG_SCHEMES, PROCESS

SCALE_BASE = os.environ.get("SCALE_BASE", f"../output/{CHANNEL}/scale_variation")


def get_results(frag, data_dir, mu_r_factor=1.0):
    """{y: [(pt, cross_section), ...]} for PROCESS and one fragmentation function."""
    if PROCESS == "sum":
        return _summed_results(data_dir, frag, mu_r_factor=mu_r_factor)
    return load_results(PROCESS, frag, data_dir=data_dir, mu_r_factor=mu_r_factor)


def scale_band(frag):
    """Minimum and maximum of the bin averages over the 7 (mu_F, mu_R) points. {} if there is no data."""
    dirs = {mu_f: CENTRAL_DIR if tag is None else f"{SCALE_BASE}/{frag}/factor_{tag}"
            for mu_f, tag in SCALE_FACTOR_TAGS.items()}
    if not all(os.path.isdir(d) for d in dirs.values()):
        return {}
    all_values = [bin_averages(get_results(frag, dirs[mu_f], mu_r)) for mu_f, mu_r in _scale_combos()]
    keys = set(all_values[0]).intersection(*all_values[1:])
    return {key: (min(v[key] for v in all_values), max(v[key] for v in all_values)) for key in keys}


def main():
    averages = {frag: bin_averages(get_results(frag, CENTRAL_DIR)) for frag, _, _, _ in FRAG_SCHEMES}
    bands = {frag: scale_band(frag) for frag, _, _, _ in FRAG_SCHEMES}
    if not any(averages.values()):
        sys.exit(f"No data found in {CENTRAL_DIR}.")

    suffix = "" if PROCESS == "sum" else f"_{PROCESS}"
    draw_panels(averages, bands, f"../plots/D0_bins_y_{CHANNEL}_{NUCLEUS}{NUCLEUS}{suffix}.pdf", ylabel_x=0.011,
                left=0.035, legend_side="split")   # the tick numbers are wider here: more room for the axis title


if __name__ == "__main__":
    main()
