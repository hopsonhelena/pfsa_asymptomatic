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

# Set input and output paths
BAM_LIST="${OUTPUT_DIR}/infected_bams.txt"

mkdir -p "$PARASITE_DIR"

# Process each parasite amplicon
for amp in "${PF_AMPS[@]}"; do
    bed="${BED_DIR}/${amp}.bed"
    raw="${PARASITE_DIR}/${amp}.pf.raw.vcf.gz"
    sorted="${PARASITE_DIR}/${amp}.pf.sort.vcf.gz"
    filtered="${PARASITE_DIR}/${amp}.pf.filtered.vcf.gz"
    dp10="${PARASITE_DIR}/${amp}.pf.filtered.dp10.vcf.gz"
    names="${PARASITE_DIR}/${amp}.sample_names.txt"

    echo "Calling parasite variants: $amp"

    # Step 1: jointly call variants across all samples
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

    # Step 2: sort the VCF and calculate INFO tags
    bcftools sort "$raw" -Ou \
    | bcftools +fill-tags -Oz -o "$sorted" -- -t all

    # Step 3: filter variants by quality, missingness, and allele frequency
    bcftools view -i 'QUAL>20 && F_MISSING<0.25 && INFO/MAF>0.05' \
        "$sorted" -Oz -o "$filtered"

    # Step 4: set genotypes with depth below 10 to missing
    bcftools +setGT "$filtered" -Oz -o "$dp10" -- -t q -n . -i 'FMT/DP<10'

    # Step 5: shorten sample names and index
    bcftools query -l "$dp10" | sed -E 's|.*/||; s|_.*||' > "$names"
    bcftools reheader -s "$names" -o "${dp10}.tmp" "$dp10"
    mv "${dp10}.tmp" "$dp10"
    bcftools index -f "$dp10"
done

echo "Depth-filtered parasite VCFs written to: $PARASITE_DIR"

# Step 6: resolve mixed genotype calls (parasite only)
# The parasite is haploid -- a het call (0/1) means multiple strains are
# present. Require at least 2 reads per allele to count as a true mixed
# infection; otherwise resolve to the better-supported homozygous call.
for vcf in "${PARASITE_DIR}"/*.pf.filtered.dp10.vcf.gz; do

    tmp="${vcf%.vcf.gz}.mixedalt.vcf.gz"
    out="${vcf%.vcf.gz}.mixedaltref.vcf.gz"

    echo "Resolving mixed genotypes: $(basename "$vcf")"

    # het with < 2 ALT reads -> REF (0/0)
    bcftools +setGT "$vcf" -Oz -o "$tmp" -- -t q -n 0 -i 'GT="0/1" & FORMAT/AD[*:1]<2'
    # het with < 2 REF reads -> ALT (1/1)
    bcftools +setGT "$tmp" -Oz -o "$out" -- -t q -n m -i 'GT="0/1" & FORMAT/AD[*:0]<2'

    bcftools index -f "$out"
    rm -f "$tmp"

done

echo "Mixed-genotype-resolved parasite VCFs written to: $PARASITE_DIR"

# Step 7: write the genotype and allele-depth table from the
# mixed-genotype-resolved VCF, so the exported
# table reflects the resolved het calls from Step 6.
for vcf in "${PARASITE_DIR}"/*.pf.filtered.dp10.mixedaltref.vcf.gz; do

    table="${vcf%.vcf.gz}.gts.txt"

    echo "Writing genotype table: $(basename "$vcf")"

    {
        printf 'CHROM\tPOS\tREF\tALT'
        bcftools query -l "$vcf" \
        | awk '{printf "\t%s_GT\t%s_AD", $1, $1} END {print ""}'
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%AD]\n' "$vcf"
    } > "$table"

done

echo "Parasite genotype tables written to: $PARASITE_DIR"
