#!/usr/bin/env bash
# Submits steps 0-6 and outputs aligned bams and per sample/amplicon coverage summary
set -euo pipefail
cd "$(dirname "$0")"        # run from amplicon-seq/ so config.sh and logs/ resolve
mkdir -p logs

source ./config.sh          # sets SLURM_ACCOUNT, SLURM_PARTITION, N_RUNS

export SBATCH_KILL_INV_DEP=yes

SLURM_OPTS=(--account="${SLURM_ACCOUNT}" --partition="${SLURM_PARTITION}") # use when no slurm array
ARRAY_OPTS=("${SLURM_OPTS[@]}" --array="0-$((N_RUNS - 1))") # use for slurm arrays

index=$(sbatch --parsable "${SLURM_OPTS[@]}" scripts/00_index.sh)   # skips indexes that already exist
trim=$(sbatch --parsable "${ARRAY_OPTS[@]}" scripts/01_trim.sh)
fqc=$(sbatch --parsable "${ARRAY_OPTS[@]}" --dependency=afterok:$trim scripts/02_fastqc.sh)
sbatch "${SLURM_OPTS[@]}" --dependency=afterok:$fqc scripts/02b_multiqc.sh
aln=$(sbatch --parsable "${ARRAY_OPTS[@]}" --dependency=afterok:$trim:$index scripts/03_align.sh)
pileup=$(sbatch --parsable "${ARRAY_OPTS[@]}" --dependency=afterok:$aln scripts/04_pileup.sh)
meancov=$(sbatch --parsable "${ARRAY_OPTS[@]}" --dependency=afterok:$pileup scripts/05_meancoverage.sh)
covtable=$(sbatch --parsable "${SLURM_OPTS[@]}" --dependency=afterok:$meancov scripts/06_coveragetable.sh)

echo "Submitted: index=$index trim=$trim fastqc=$fqc align=$aln pileup=$pileup meancov=$meancov covtable=$covtable"
echo "When $covtable finishes, run steps 07-10 manually to select cutoffs (see README)."
