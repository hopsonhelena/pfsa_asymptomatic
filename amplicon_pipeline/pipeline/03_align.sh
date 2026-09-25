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
BAM_DIR="${WORK_DIR}/${run}/${BAM_SUBDIR}"

mkdir -p "$BAM_DIR"

echo "Aligning run: $run"

for r2 in "${TRIMMED_DIR}/${run%R}X"*_R2_001.cutadapt.fastq.gz; do

    sample=$(basename "$r2" _R2_001.cutadapt.fastq.gz)

    r1="${TRIMMED_DIR}/${sample}_R1_001.cutadapt.fastq.gz"
    bam="${BAM_DIR}/${sample}.joint.human.pf.sorted.bam"

    echo "Aligning sample: $sample"

    bwa mem -M -t "$SLURM_CPUS_PER_TASK" "$REF" \
        "$r1" \
        "$r2" \
        | samtools sort -o "$bam"

    samtools index "$bam"

done

echo "Alignment complete: $run"
