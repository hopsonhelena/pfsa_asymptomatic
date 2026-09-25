#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --out=logs/s-%x.%j.out
#SBATCH --error=logs/s-%x.%j.err

set -euo pipefail

source "./config.sh"
set +u
eval "$(conda shell.bash hook)"
conda activate amplicon-pipeline
set -u

MULTIQC_DIR="${OUTPUT_DIR}/multiqc"
mkdir -p "$MULTIQC_DIR"

fastqc_dirs=()

for run in "${ALL_RUNS[@]}"; do
    fastqc_dir="${WORK_DIR}/${run}/${FASTQC_SUBDIR}"

    if [[ ! -d "$fastqc_dir" ]]; then
        echo "FastQC directory not found for run: $run" >&2
        exit 1
    fi

    fastqc_dirs+=("$fastqc_dir")
done


echo "Creating MultiQC report across all runs"

multiqc "${fastqc_dirs[@]}" \
    --outdir "$MULTIQC_DIR" \
    --force

echo "MultiQC report written: ${MULTIQC_DIR}/multiqc_report.html"

