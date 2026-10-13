#!/usr/bin/env bash

# Builds bam lists used for variant calling 

set -euo pipefail 

source "./config.sh"

INFECTION_CALLS="${OUTPUT_DIR}/infection_calling/infection_calls.csv"
AVAILABLE_BAMS="${OUTPUT_DIR}/available_bams.txt"
INFECTED_BAMS="${OUTPUT_DIR}/infected_bams.txt"
ALL_BAMS="${OUTPUT_DIR}/all_bams.txt"


if [ ! -f "$INFECTION_CALLS" ]; then
    echo "ERROR: $INFECTION_CALLS not found. Run step 10 first." >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
> "$AVAILABLE_BAMS"

for run in "${ALL_RUNS[@]}"; do
    bam_dir="${WORK_DIR}/${run}/${BAM_SUBDIR}"
    find "$bam_dir" -name '*.bam' | sort >> "$AVAILABLE_BAMS"
done


# Match "/<sample>_" against the BAM paths. The leading "/" and trailing "_"
# stop 24118X10 from also matching 24118X100.
awk -F, 'NR > 1 {print "/" $1 "_"}' "$INFECTION_CALLS" |
    grep -Ff - "$AVAILABLE_BAMS" > "$ALL_BAMS" || true

awk -F, 'NR > 1 && $NF == "INFECTED" {print "/" $1 "_"}' "$INFECTION_CALLS" |
    grep -Ff - "$AVAILABLE_BAMS" > "$INFECTED_BAMS" || true

rm -f "$AVAILABLE_BAMS"

if [ ! -s "$ALL_BAMS" ]; then
    echo "ERROR: no BAMs matched the samples in $INFECTION_CALLS" >&2
    exit 1
fi

n_samples=$(( $(wc -l < "$INFECTION_CALLS") - 1 ))
echo "Samples in infection calls: $n_samples"
echo "All BAMs:      $(wc -l < "$ALL_BAMS")  -> $ALL_BAMS"
echo "Infected BAMs: $(wc -l < "$INFECTED_BAMS")  -> $INFECTED_BAMS"

