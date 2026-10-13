#!/usr/bin/env bash
# Runs steps 0-6 in order locally. Use instead of submit_alignment.sh
set -euo pipefail
cd "$(dirname "$0")"        # run from amplicon-seq/ so config.sh resolves

source ./config.sh          # sets ALL_RUNS

bash scripts/00_index.sh

for i in "${!ALL_RUNS[@]}"; do
    for step in 01_trim 02_fastqc 03_align 04_pileup 05_meancoverage; do
        bash "scripts/${step}.sh" "$i"
    done
done

bash scripts/02b_multiqc.sh
bash scripts/06_coveragetable.sh

echo "Next, run steps 07-10 manually to select cutoffs (see README)."
