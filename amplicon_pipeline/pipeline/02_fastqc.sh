#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
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

TRIMMED_DIR="${WORK_DIR}/${run}/${TRIMMED_SUBDIR}"
FASTQC_DIR="${WORK_DIR}/${run}/${FASTQC_SUBDIR}"

mkdir -p "$FASTQC_DIR" 

echo "FastQC run: $run"

fastqc -t "$SLURM_CPUS_PER_TASK" \
    -o "$FASTQC_DIR" \
    "${TRIMMED_DIR}"/*_R1_001.cutadapt.fastq.gz \
    "${TRIMMED_DIR}"/*_R2_001.cutadapt.fastq.gz

echo "FastQC complete: $run" 
