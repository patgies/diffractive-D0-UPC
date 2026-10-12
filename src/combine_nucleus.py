import os
import sys

# Writes tables of the final D0 cross sections dsigma/(dy dpD0) [mb/GeV]: the run_nucleus.sh output
# integrated over the impact parameter b_d, with the prefactors, as in the plots (plotting_scripts/D0.py).
# One file per process and rapidity, in output/<CHANNEL>/combined/.
# usage (from src/ or plotting_scripts/): CHANNEL=An0n FRAG_TYPE=HymnD python3 ../src/combine_nucleus.py

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "plotting_scripts"))
from D0 import load_results, NUCLEUS, CHANNEL, CENTRAL_DIR, OUTPUT_BASE, FLUX_TAG

FRAG = os.environ.get("FRAG_TYPE", "HymnD")
COMBINED_DIR = os.environ.get("COMBINED_DIR", f"{OUTPUT_BASE}/{CHANNEL}/combined")


def write_combined(process, results):
    os.makedirs(COMBINED_DIR, exist_ok=True)
    for y, points in sorted(results.items()):
        ytag = f"{y:.1f}".replace(".", "")
        outfile = f"{COMBINED_DIR}/combined_d0_{process}_{FRAG}_{CHANNEL}_{NUCLEUS}{FLUX_TAG}_y{ytag}.dat"
        with open(outfile, "w") as f:
            f.write(f"# {process} D0 cross section, {NUCLEUS} target, {FRAG} fragmentation, "
                    f"{CHANNEL} channel, integrated over the impact parameter b_d\n")
            f.write(f"# fixed rapidity y : {y}\n")
            f.write("# pD0  dsigma_dy_dpD0 [mb/GeV]\n")
            for pt, cross_section in points:
                f.write(f"{pt:.6g}  {cross_section:.6g}\n")
        print(f"Wrote {outfile}")


def main():
    results = {process: load_results(process, FRAG) for process in ("exclusive", "diffractive")}
    if not any(results.values()):
        sys.exit(f"No files found for NUCLEUS={NUCLEUS}, FRAG_TYPE={FRAG}, CHANNEL={CHANNEL} "
                 f"in {CENTRAL_DIR} -- run ../run_scripts/run_nucleus.sh first.")
    for process, process_results in results.items():
        write_combined(process, process_results)


if __name__ == "__main__":
    main()
