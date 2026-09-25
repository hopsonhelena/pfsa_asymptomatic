#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/s-%x.%j.out
#SBATCH --error=logs/s-%x.%j.err

set -euo pipefail

source "./config.sh"
set +u
eval "$(conda shell.bash hook)"
conda activate amplicon-pipeline
set -u

BAM_LIST="${OUTPUT_DIR}/all_bams.txt"
mkdir -p "$HUMAN_DIR"

# Process each human amplicon
for amp in "${HUMAN_AMPS[@]}"; do
    bed="${BED_DIR}/${amp}.bed"
    raw="${HUMAN_DIR}/${amp}.human.raw.vcf.gz"
    sorted="${HUMAN_DIR}/${amp}.human.sort.vcf.gz"
    filtered="${HUMAN_DIR}/${amp}.human.filtered.vcf.gz"
    dp10="${HUMAN_DIR}/${amp}.human.filtered.dp10.vcf.gz"
    names="${HUMAN_DIR}/${amp}.sample_names.txt"
    table="${HUMAN_DIR}/${amp}.human.filtered.dp10.gts.txt"

    echo "Calling human variants: $amp"

    # Jointly call variants across all samples
    bcftools mpileup \
        --max-depth 1000000 \
        --max-idepth 100000 \
        -Ou \
        -f "$REF" \
        -b "$BAM_LIST" \
        --per-sample-mF \
        --regions-file "$bed" \
        --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    | bcftools call -m -v -Oz -o "$raw"

    # Sort the VCF and calculate INFO tags
    bcftools sort "$raw" -Ou \
    | bcftools +fill-tags -Oz -o "$sorted" -- -t all

    # Filter variants by quality, missingness, and allele frequency
    bcftools view -i 'QUAL>20 && F_MISSING<0.25 && INFO/MAF>0.05' \
        "$sorted" -Oz -o "$filtered"

    # Set genotypes with depth below 10 to missing
    bcftools +setGT "$filtered" -Oz -o "$dp10" -- -t q -n . -i 'FMT/DP<10'

    # Shorten sample names and index the final VCF
    bcftools query -l "$dp10" | sed -E 's|.*/||; s|_.*||' > "$names"
    bcftools reheader -s "$names" -o "${dp10}.tmp" "$dp10"
    mv "${dp10}.tmp" "$dp10"
    bcftools index -f "$dp10"

    # Write the genotype and allele-depth table
    {
        printf 'CHROM\tPOS\tREF\tALT'
        bcftools query -l "$dp10" \
        | awk '{printf "\t%s_GT\t%s_AD", $1, $1} END {print ""}'
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%AD]\n' "$dp10"
    } > "$table"
done

echo "Human VCFs and genotype tables written to: $HUMAN_DIR"

