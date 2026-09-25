"""
Assess plate contamination using negative and positive control coverage.

For each amplicon on each plate, compute:

    control_ratio = average negative-control coverage / average positive-control coverage

A high ratio is a sign of contamination. A plate is flagged for exclusion if
more than half of its amplicons have a ratio above some cutoff.

Run this script in two steps: 

  Step 1 - run with no --cutoff. This just plots the ratios so you can
  look at them and decide where a sensible cutoff is.

      python 07_plateQC.py --cov coverage.txt --meta metadata.csv

  Step 2 - once you've picked a cutoff (e.g. 0.02), rerun with it to get
  the actual PASS/FAIL and KEEP/EXCLUDE calls.

      python 07_plateQC.py --cov coverage.txt --meta metadata.csv --cutoff 0.02
"""

import os
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


# Fraction of amplicons that must pass
# A plate is excluded if more than this fraction of its amplicons have 
# coverage ratio over empirical threshold 
MAX_FAILED_FRACTION = 0.50


# Input and output file paths, if --cov/--meta/--results-dir is unspecified
OUTPUT_DIR = Path(os.environ['OUTPUT_DIR'])
DEFAULT_COVERAGE = Path(os.environ['COVERAGE_TABLE'])
DEFAULT_METADATA = Path(os.environ['METADATA'])
DEFAULT_RESULTS_DIR = OUTPUT_DIR / 'plateQC'


# Step 1: compute the control_ratio for every plate/amplicon combination
def compute_control_ratios(coverage_df, metadata_df, cutoff=None):
    """
    Combine the coverage table and the metadata table, then compute the
    negative/positive control ratio for every (plate, amplicon) pair.

    cutoff: if given, also mark each row as PASS or FAIL
    """

    # Clean up the coverage table 
    coverage_df = coverage_df.rename(columns={"sample": "seq_sample_id"})
    coverage_df["seq_sample_id"] = coverage_df["seq_sample_id"].astype(str)

    # Clean up the metadata table
    metadata_df = metadata_df.copy()
    metadata_df["seq_sample_id"] = metadata_df["seq_sample_id"].astype(str)
    metadata_df["plateID"] = metadata_df["plateID"].astype(str)
    metadata_df["sample_type"] = metadata_df["sample_type"].str.lower()
    metadata_df["control_type"] = metadata_df["control_type"].str.lower()

    # Keep only the control samples 
    controls_df = metadata_df[metadata_df["sample_type"] == "control"]

    # merge control metadata with their coverage values 
    merged_df = controls_df.merge(coverage_df, on="seq_sample_id", how="inner")

    # For each plate + amplicon, compute the average coverage of
    # the negative controls and the average coverage of the
    # positive controls

    results = []  # list of dictionaries, one per plate, amplicon

    # group by plate and amplicon
    grouped = merged_df.groupby(["plateID", "amplicon"])

    for (plate_id, amplicon), group in grouped:
        negative_rows = group[group["control_type"] == "negative"]
        positive_rows = group[group["control_type"] == "positive"]

        negative_mean = negative_rows["mean_coverage"].mean()
        positive_mean = positive_rows["mean_coverage"].mean()

        # Guard against dividing by zero or missing controls
        if pd.isna(negative_mean) or pd.isna(positive_mean) or positive_mean == 0:
            ratio = None
        else:
            ratio = negative_mean / positive_mean

        row = {
            "plateID": plate_id,
            "amplicon": amplicon,
            "negative": negative_mean,
            "positive": positive_mean,
            "control_ratio": ratio,
        }

        if cutoff is not None:
            if ratio is not None and ratio <= cutoff:
                row["amplicon_status"] = "PASS"
            else:
                row["amplicon_status"] = "FAIL"

        results.append(row)

    ratios_df = pd.DataFrame(results)
    ratios_df = ratios_df.sort_values(["plateID", "amplicon"]).reset_index(drop=True)
    return ratios_df


# Step 2: combine per-amplicon PASS/FAIL calls to a per-plate decision

def summarize_plates(ratios_df):
    """
    For each plate, count how many amplicons passed vs failed, and decide
    whether the whole plate should be kept or excluded.
    """

    summary_rows = []

    for plate_id, group in ratios_df.groupby("plateID"):
        n_amplicons = len(group)
        n_failed = (group["amplicon_status"] == "FAIL").sum()
        frac_failed = n_failed / n_amplicons

        if frac_failed > MAX_FAILED_FRACTION:
            plate_status = "EXCLUDE"
            reason = f"{n_failed}/{n_amplicons} amplicons failed control-ratio QC"
        else:
            plate_status = "KEEP"
            reason = "Passed control-ratio QC"

        summary_rows.append({
            "plateID": plate_id,
            "n_amplicons": n_amplicons,
            "n_failed": n_failed,
            "frac_failed": frac_failed,
            "plate_status": plate_status,
            "reason": reason,
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df = summary_df.sort_values("plateID").reset_index(drop=True)
    return summary_df


# Plotting

def plot_control_ratios(ratios_df, output_path, cutoff=None, summary_df=None):
    """
    Scatter plot of control_ratio for every amplicon, grouped by plate.
    If given, cutoff represented by dashed line.
    """

    plates = sorted(ratios_df["plateID"].unique())
    # give each plate a number 
    plate_to_x = {plate: i for i, plate in enumerate(plates)}

    fig, ax = plt.subplots(figsize=(10, 6))

    # plot each amplicon as its own set of points, in its own color
    for amplicon, group in ratios_df.groupby("amplicon"):
        x_values = group["plateID"].map(plate_to_x)
        y_values = group["control_ratio"]
        ax.scatter(x_values, y_values, label=amplicon, alpha=0.8)

    # shade excluded plates
    if summary_df is not None:
        excluded_plates = summary_df[summary_df["plate_status"] == "EXCLUDE"]["plateID"]
        for plate in excluded_plates:
            x = plate_to_x[plate]
            ax.axvspan(x - 0.45, x + 0.45, color="red", alpha=0.08)

    # draw the cutoff line
    if cutoff is not None:
        ax.axhline(cutoff, color="red", linestyle="--", label="cutoff")

    ax.set_yscale("log")
    ax.set_xticks(range(len(plates)))
    ax.set_xticklabels(plates, rotation=45, ha="right")
    ax.set_xlabel("Plate")
    ax.set_ylabel("Negative / positive control coverage ratio")
    ax.set_title("Control contamination QC by plate")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_ratio_elbow(ratios_df, output_path):
    """
    Sorts every control_ratio from smallest to largest and plots it. 
    """

    # drop rows where we couldn't compute a ratio
    valid_ratios = ratios_df.dropna(subset=["control_ratio"])
    sorted_ratios = valid_ratios.sort_values("control_ratio").reset_index(drop=True)

    # x-axis is just 1, 2, 3, ... one number per data point
    x_values = range(1, len(sorted_ratios) + 1)
    y_values = sorted_ratios["control_ratio"]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(x_values, y_values, marker="o", markersize=3)

    # reference lines at 0.00, 0.01, 0.02, 0.03, 0.04, 0.05
    for reference_value in [0.00, 0.01, 0.02, 0.03, 0.04, 0.05]:
        ax.axhline(reference_value, color="gray", linewidth=0.5, linestyle="--", zorder=0)

    ax.set_xlabel("Rank (amplicon x plate, sorted by ratio)")
    ax.set_ylabel("Negative / positive control coverage ratio")
    ax.set_title("Control ratio elbow plot")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


# Main program

if __name__ == "__main__":

    #  command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--cov", type=Path, default=DEFAULT_COVERAGE,
                         help="Path to the coverage file")
    parser.add_argument("--meta", type=Path, default=DEFAULT_METADATA,
                         help="Path to the sample metadata CSV")
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR,
                         help="Folder to write output files into")
    parser.add_argument("--cutoff", type=float, default=None,
                         help="Ratio cutoff for PASS/FAIL. Leave blank to just "
                              "make plots and pick a cutoff by eye.")
    args = parser.parse_args()

    # load the input files 
    coverage_df = pd.read_csv(args.cov, sep="\t")
    metadata_df = pd.read_csv(args.meta)

    # control ratios computation
    ratios_df = compute_control_ratios(coverage_df, metadata_df, cutoff=args.cutoff)

    # save the ratio table, sorted lowest ratio to highest 
    args.results_dir.mkdir(parents=True, exist_ok=True)
    ratios_path = args.results_dir / "control_coverage_ratios.csv"
    ratios_df_sorted = ratios_df.sort_values("control_ratio")
    ratios_df_sorted.to_csv(ratios_path, index=False)

    scatter_plot_path = args.results_dir / "control_coverage_ratios.png"

    if args.cutoff is None:

        plot_control_ratios(ratios_df, scatter_plot_path)

        elbow_plot_path = args.results_dir / "control_ratio_elbow.png"
        plot_ratio_elbow(ratios_df, elbow_plot_path)

        print(f"Ratio table written to: {ratios_path}")
        print(f"Scatter plot written to: {scatter_plot_path}")
        print(f"Elbow plot written to: {elbow_plot_path}")
        print("No --cutoff given. Look at the plots, then rerun with --cutoff.")

    else:
        # with cutoff, compute plate-level decisions.
        summary_df = summarize_plates(ratios_df)

        summary_path = args.results_dir / "plate_qc_summary.csv"
        summary_df.to_csv(summary_path, index=False)

        plot_control_ratios(ratios_df, scatter_plot_path, cutoff=args.cutoff, summary_df=summary_df)

        print(summary_df[["plateID", "plate_status", "reason"]].to_string(index=False))
        print()
        print(f"Ratio table written to: {ratios_path}")
        print(f"Plate summary written to: {summary_path}")
        print(f"Scatter plot written to: {scatter_plot_path}")
