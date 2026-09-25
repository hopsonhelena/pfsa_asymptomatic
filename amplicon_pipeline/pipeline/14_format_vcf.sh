#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/s-%x.%j.out
#SBATCH --error=logs/s-%x.%j.err

set -euo pipefail

source "./config.sh"
set +u
eval "$(conda shell.bash hook)"
conda activate amplicon-pipeline
set -u

# Genotype/allele-depth table
write_table() {
    {
        printf 'CHROM\tPOS\tREF\tALT'
        bcftools query -l "$1" | awk '{printf "\t%s_GT\t%s_AD", $1, $1} END {print ""}'
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%AD]\n' "$1"
    } > "$2"
}

# List of Febrile samples to exclude (optional): read the seq_sample_id column of FEBRILE_LIST, one ID per line.
FEBRILE_IDS="${OUTPUT_DIR}/febrile_sample_ids.txt"
: > "$FEBRILE_IDS"
if [[ -f "${FEBRILE_LIST:-}" ]]; then
    tr -d '"\r' < "$FEBRILE_LIST" \
    | awk -F, 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "seq_sample_id") col = i; next}
               col {print $col}
               END {if (!col) {print "ERROR: no seq_sample_id column in FEBRILE_LIST" > "/dev/stderr"; exit 1}}' \
    > "$FEBRILE_IDS"
fi

# Human: remove febrile samples, recompute INFO tags
for amp in "${HUMAN_AMPS[@]}"; do
    out="${HUMAN_DIR}/${amp}.human.filtered.final.vcf.gz"
    bcftools view --force-samples -S ^"$FEBRILE_IDS" \
        "${HUMAN_DIR}/${amp}.human.filtered.dp10.vcf.gz" \
    | bcftools +fill-tags -Oz -o "$out" -- -t all
    write_table "$out" "${HUMAN_DIR}/${amp}.human.filtered.final.gts.txt"
done

# Parasite: concatenate all amplicons, remove febrile samples, recompute INFO tags
PF_ALL="${PARASITE_DIR}/pf.all.filtered.final.vcf.gz"
bcftools concat -a "${PARASITE_DIR}"/*.pf.filtered.dp10.mixedaltref.vcf.gz \
| bcftools sort \
| bcftools view --force-samples -S ^"$FEBRILE_IDS" \
| bcftools +fill-tags -Oz -o "$PF_ALL" -- -t all
write_table "$PF_ALL" "${PARASITE_DIR}/pf.all.filtered.final.gts.txt"

echo "Done: $HUMAN_DIR and $PARASITE_DIR"
