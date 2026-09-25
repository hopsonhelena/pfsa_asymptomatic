#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --out=logs/s-%x.%A_%a.out
#SBATCH --error=logs/s-%x.%A_%a.err

set -euo pipefail

source "./config.sh"
set +u
eval "$(conda shell.bash hook)"
conda activate amplicon-pipeline
set -u

# Select run from slurm arrray
run=${ALL_RUNS[$SLURM_ARRAY_TASK_ID]}

FASTQ_DIR="${DATA_DIR}/${run}/${FASTQ_SUBDIR}"
TRIMMED_DIR="${WORK_DIR}/${run}/${TRIMMED_SUBDIR}"

mkdir -p "$TRIMMED_DIR"

echo "Trimming run: $run"

for r1 in "${FASTQ_DIR}/${run%R}X"*_R1_001.fastq.gz; do

    sample=$(basename "$r1" _R1_001.fastq.gz)
    r2="${FASTQ_DIR}/${sample}_R2_001.fastq.gz"

    cutadapt \
        -j "$SLURM_CPUS_PER_TASK" \
        -m "$MIN_LENGTH" \
        -a "$ADAPTER_R1" \
        -A "$ADAPTER_R2" \
        -o "${TRIMMED_DIR}/${sample}_R1_001.cutadapt.fastq.gz" \
        -p "${TRIMMED_DIR}/${sample}_R2_001.cutadapt.fastq.gz" \
        "$r1" "$r2"

done

echo "Trim complete: $run"
