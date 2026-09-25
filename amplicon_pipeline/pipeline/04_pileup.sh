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

BAM_DIR="${WORK_DIR}/${run}/${BAM_SUBDIR}"
DEPTH_DIR="${WORK_DIR}/${run}/${DEPTH_SUBDIR}"
BED_FILE="${BED_DIR}/human_pf_amplicons.bed"

mkdir -p "$DEPTH_DIR"

echo "Creating pileups for run: $run"

for bam in "${BAM_DIR}/${run%R}X"*.joint.human.pf.sorted.bam; do

    # Keep only sample name (e.g. 24118X100) from the full bam name
    sample=$(basename "$bam" | cut -d '_' -f1)
    pileup="${DEPTH_DIR}/${sample}.pileup.txt"

    echo "Creating pileup for sample: $sample"

    samtools mpileup \
        -f "$REF" \
        -l "$BED_FILE" \
        -aa \
        -Q 1 \
        -o "$pileup" \
        "$bam"

done

echo "Pileup complete: $run"



