#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --out=logs/s-%x.%A_%a.out
#SBATCH --error=logs/s-%x.%A_%a.err

set -euo pipefail

source "./config.sh"
set +u
eval "$(conda shell.bash hook)"
conda activate amplicon-pipeline
set -u

# Select run from slurm array
run=${ALL_RUNS[$SLURM_ARRAY_TASK_ID]}

DEPTH_DIR="${WORK_DIR}/${run}/${DEPTH_SUBDIR}"
BED_FILE="${BED_DIR}/human_pf_amplicons.bed"

echo "Calculating mean coverage for run: $run"

for pileup in "${DEPTH_DIR}"/*.pileup.txt; do

    sample=$(basename "$pileup" .pileup.txt)
    outfile="${DEPTH_DIR}/${sample}.mean_coverage.txt"
    > "$outfile"

    echo "Calculating mean coverage for sample: $sample"

    while IFS=$'\t' read -r chrom start end amplicon rest; do
        mean=$(awk -v c="$chrom" -v s="$start" -v e="$end" \
            '$1==c && $2>=s && $2<=e {sum+=$4; n++} END {print (n>0 ? sum/n : 0)}' \
            "$pileup")
        printf '%s\t%s\t%s\n' "$amplicon" "$sample" "$mean" >> "$outfile"
    done < "$BED_FILE"

done

echo "Mean coverage complete: $run"

