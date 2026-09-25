#!/usr/bin/env bash
# Submits steps 11-14 (variant calling). Run after steps 07-10.
set -euo pipefail
cd "$(dirname "$0")"        # run from pipeline directory so config.sh and logs/ resolve
mkdir -p logs

source ./config.sh          # sets SLURM_ACCOUNT, SLURM_PARTITION, N_RUNS

export SBATCH_KILL_INV_DEP=yes

SLURM_OPTS=(--account="${SLURM_ACCOUNT}" --partition="${SLURM_PARTITION}") # use when no slurm array


bash pipeline/11_list_bams.sh
human=$(sbatch "${SLURM_OPTS[@]}" --parsable pipeline/12_call_human_variants.sh)
parasite=$(sbatch "${SLURM_OPTS[@]}" --parsable pipeline/13_call_parasite_variants.sh)
sbatch "${SLURM_OPTS[@]}" --dependency=afterok:$human:$parasite pipeline/14_format_vcf.sh

echo "Submitted: human=$human parasite=$parasite"
