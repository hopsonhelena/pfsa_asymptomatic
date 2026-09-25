#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --out=logs/s-%x.%j.out
#SBATCH --error=logs/s-%x.%j.err

# Output table summarizing coverage across all amplicons and samples 

set -euo pipefail
source "./config.sh"
set +u
eval "$(conda shell.bash hook)"
conda activate amplicon-pipeline
set -u

mkdir -p $OUTPUT_DIR

echo "Combining mean coverage files across all runs"

printf 'amplicon\tsample\tmean_coverage\n' > "$COVERAGE_TABLE"

for run in "${ALL_RUNS[@]}"; do

    DEPTH_DIR="${WORK_DIR}/${run}/${DEPTH_SUBDIR}"
    coverage_files=("${DEPTH_DIR}"/*.mean_coverage.txt)

    echo "Adding mean coverage files for run: $run"

    for coverage_file in "${coverage_files[@]}"; do
        cat "$coverage_file" >> "$COVERAGE_TABLE"
    done

done

total_rows=$(($(wc -l < "$COVERAGE_TABLE") - 1))

echo "Coverage table written: $COVERAGE_TABLE"
echo "Total rows (excluding header): $total_rows"

