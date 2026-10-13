#!/usr/bin/env bash
# Runs steps 11-14 (variant calling) in order locally. Use instead of submit_variant_calling.sh
set -euo pipefail
cd "$(dirname "$0")"        # run from amplicon-seq/ so config.sh resolves

bash scripts/11_list_bams.sh
bash scripts/12_call_human_variants.sh
bash scripts/13_call_parasite_variants.sh
bash scripts/14_format_vcf.sh
