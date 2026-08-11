#!/usr/bin/env bash
#
# run_snptest_1701_infected_samples.sh
#
# Test for association between HbS genotype and phenotypes in infected samples: 
# complexity of infection (continuous and binary) and parasitemia 
#
# Requirements: SNPTEST v2.5.2 (https://www.chg.ox.ac.uk/~gav/snptest/)
# set SNPTEST_BIN below to its path.
 
set -euo pipefail

## CONFIG
SNPTEST_BIN="snptest_v2.5.2"
SAMPLES='cameroon_1701_infected_samples.sample'
VCF='cameroon_hbb_genotypes_1701_infected_samples.vcf.gz'
OUTDIR="results/snptest"

mkdir -p "$OUTDIR"/logs

# run_snptest PHENOTYPE [COVARIATE...]
# Runs one SNPTEST association test. Any covariates given are passed
# straight to -cov_names; omit them for the unadjusted run.
run_snptest () {
    pheno="$1"
    shift
    covariates="$@"
 
    if [ -n "$covariates" ]; then
        tag="${pheno}_$(echo $covariates | tr ' ' '_')"
        "$SNPTEST_BIN" -data "$VCF" "$SAMPLES" -genotype_field GT -pheno "$pheno" \
            -frequentist 1 -method threshold -cov_names $covariates \
            -o "$OUTDIR/$tag.out" > "$OUTDIR/logs/$tag.log" 2>&1
    else
        tag="$pheno"
        "$SNPTEST_BIN" -data "$VCF" "$SAMPLES" -genotype_field GT -pheno "$pheno" \
            -frequentist 1 -method threshold \
            -o "$OUTDIR/$tag.out" > "$OUTDIR/logs/$tag.log" 2>&1
    fi
}


## ANALYSES RUN 
# Complexity of infection (continuous), no covariates
run_snptest AMA1_coi
run_snptest SERA2_coi
 
# Complexity of infection (continuous), adjusted for village + season
run_snptest AMA1_coi  Village Season
run_snptest SERA2_coi Village Season
 
# Complexity of infection (binary), no covariates 
run_snptest AMA1_poly
run_snptest SERA2_poly

# Complexity of infection (binary), adjusted for village + season
run_snptest AMA1_poly  Village Season
run_snptest SERA2_poly Village Season
 
# Parasitemia (continuous), no covariates 
run_snptest Parasitemia_0_4

# Parasitemia (continuous), adjusted for village + season
run_snptest Parasitemia_0_4 Village Season 
