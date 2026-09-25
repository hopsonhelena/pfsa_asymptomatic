#!/usr/bin/env bash
# Submits steps 1-6 and outputs aligned bams and per sample/amplicon coverage summary
set -euo pipefail
cd "$(dirname "$0")"        # run from pipeline directory so config.sh and logs/ resolve
mkdir -p logs

source ./config.sh          # sets SLURM_ACCOUNT, SLURM_PARTITION, N_RUNS

export SBATCH_KILL_INV_DEP=yes

SLURM_OPTS=(--account="${SLURM_ACCOUNT}" --partition="${SLURM_PARTITION}") # use when no slurm array
ARRAY_OPTS=("${SLURM_OPTS[@]}" --array="0-$((N_RUNS - 1))") # use for slurm arrays

trim=$(sbatch --parsable "${ARRAY_OPTS[@]}" pipeline/01_trim.sh)
fqc=$(sbatch --parsable "${ARRAY_OPTS[@]}" --dependency=afterok:$trim pipeline/02_fastqc.sh)
sbatch "${SLURM_OPTS[@]}" --dependency=afterok:$fqc pipeline/02b_multiqc.sh
aln=$(sbatch --parsable "${ARRAY_OPTS[@]}" --dependency=afterok:$trim pipeline/03_align.sh)
pileup=$(sbatch --parsable "${ARRAY_OPTS[@]}" --dependency=afterok:$aln pipeline/04_pileup.sh)
meancov=$(sbatch --parsable "${ARRAY_OPTS[@]}" --dependency=afterok:$pileup pipeline/05_meancoverage.sh)
covtable=$(sbatch --parsable "${SLURM_OPTS[@]}" --dependency=afterok:$meancov pipeline/06_coveragetable.sh)

echo "Submitted: trim=$trim fastqc=$fqc align=$aln pileup=$pileup meancov=$meancov covtable=$covtable"
echo "When $covtable finishes, run steps 07-10 manually to select cutoffs (see README)."
