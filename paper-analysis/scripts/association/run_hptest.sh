#!/usr/bin/env bash
# Test associations between human HbS genotype and P. falciparum variants,
# with and without adjusting for village and season.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname -- "$(dirname -- "$SCRIPT_DIR")")"
cd "$REPO_ROOT"

## CONFIG
HPTEST_BIN="hptest_v2.2.0" # set to hptest path 
SAMPLES="data/hptest_samples.txt"
HBB="data/cameroon_hbb_genotypes_1701_infected_samples.vcf.gz"
PF="data/cameroon_Pf_genotypes_1701_infected_samples.vcf.gz"
OUTDIR="results/association/hptest"
mkdir -p "$OUTDIR"/logs
  
# run_hptest [COVARIATE...]
# Runs one hptest association test. Any covariates given are passed
# straight to -covariates; omit them for the unadjusted run.
run_hptest () {
    local covariates=("$@")

    if (( ${#covariates[@]} )); then
        local covariate_tag
        covariate_tag="$(IFS=_; echo "${covariates[*]}")"
        tag="hbb_vs_Pf_${covariate_tag}"
        "$HPTEST_BIN" -outcome "$PF" -predictor "$HBB" -s "$SAMPLES" \
            -treat-outcome-as-haploid -output-parameters all -output-counts \
            -covariates "${covariates[@]}" \
            -o "$OUTDIR/$tag.hptest.csv" > "$OUTDIR/logs/$tag.log" 2>&1
    else
        tag="hbb_vs_Pf_no_covariates"
        "$HPTEST_BIN" -outcome "$PF" -predictor "$HBB" -s "$SAMPLES" \
            -treat-outcome-as-haploid -output-parameters all -output-counts \
            -o "$OUTDIR/$tag.hptest.csv" > "$OUTDIR/logs/$tag.log" 2>&1
    fi
}
 

## ANALYSES RUN
run_hptest
run_hptest Village Season
