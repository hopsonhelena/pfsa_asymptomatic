#!/usr/bin/env bash
#
# run_hptest.sh
#
# hptest v2.2.0 host-parasite association tests between HbS genotype
# and P. falciparum variants, with and without adjusting for Village and Season.
#
# Requirements: hptest v2.2.0 which is included as part of the qctool package (download here: https://enkre.net/code/)
 
set -euo pipefail
 
## CONFIG
HPTEST_BIN="hptest_v2.2.0"
SAMPLES="data/hptest_samples.txt"
HBB="data/cameroon_hbb_genotypes_1701_infected_samples.vcf.gz"
PF="data/cameroon_Pf_genotypes_1701_infected_samples.vcf.gz"
OUTDIR="results/hptest"
mkdir -p "$OUTDIR"/logs
  
# run_hptest [COVARIATE...]
# Runs one hptest association test. Any covariates given are passed
# straight to -covariates; omit them for the unadjusted run.
run_hptest () {
    covariates="$@"
 
    if [ -n "$covariates" ]; then
        tag="hbb_vs_Pf_$(echo $covariates | tr ' ' '_')"
        "$HPTEST_BIN" -outcome "$PF" -predictor "$HBB" -s "$SAMPLES" \
            -treat-outcome-as-haploid -output-parameters all -output-counts \
            -covariates $covariates \
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
