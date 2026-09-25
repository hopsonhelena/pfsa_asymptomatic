#!/usr/bin/env bash
# Test associations between HbS genotype and 
# P. falciparum infection status and gametocyte presence
# with and without covariates 

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname -- "$(dirname -- "$SCRIPT_DIR")")"
cd "$REPO_ROOT"

### CONFIG 
SNPTEST_BIN="snptest_v2.5.2" # set to snptest path
SAMPLES='data/cameroon_2246_samples.sample'
VCF='data/cameroon_hbb_genotypes_2246_samples.vcf.gz'
OUTDIR="results/association/snptest"

mkdir -p "$OUTDIR"/logs

# run_snptest PHENOTYPE [COVARIATE...]
# Runs one SNPTEST association test. Any covariates given are passed
# straight to -cov_names; omit them for the unadjusted run.
run_snptest () {
    pheno="$1"
    shift
    local covariates=("$@")

    if (( ${#covariates[@]} )); then
        local covariate_tag
        covariate_tag="$(IFS=_; echo "${covariates[*]}")"
        tag="${pheno}_${covariate_tag}"
        "$SNPTEST_BIN" -data "$VCF" "$SAMPLES" -genotype_field GT -pheno "$pheno" \
            -frequentist 1 -method threshold -cov_names "${covariates[@]}" \
            -o "$OUTDIR/$tag.out" > "$OUTDIR/logs/$tag.log" 2>&1
    else
        tag="$pheno"
        "$SNPTEST_BIN" -data "$VCF" "$SAMPLES" -genotype_field GT -pheno "$pheno" \
            -frequentist 1 -method threshold \
            -o "$OUTDIR/$tag.out" > "$OUTDIR/logs/$tag.log" 2>&1
    fi
}


## ANALYSES RUN
# Infection status from sequencing (binary), no covariates
run_snptest Infected

# Infection status from sequencing (binary), adjusted for village + season
run_snptest Infected Village Season

# Infection status from microscopy (binary), no covariates
run_snptest Infected_microscopy

# Infection status from microscopy (binary), adjusted for village + season
run_snptest Infected_microscopy Village Season
 
# Gametocyte presence from microscopy (binary), no covariates 
run_snptest Gametocytes_YN

# Gametocyte presence from microscopy (binary), adjusted for village + season
run_snptest Gametocytes_YN Village Season 
