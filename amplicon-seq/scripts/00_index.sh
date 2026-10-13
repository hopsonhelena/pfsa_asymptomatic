#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --out=logs/s-%x.%j.out
#SBATCH --error=logs/s-%x.%j.err

# Create the BWA and samtools indexes for the reference, if they do not already exist

set -euo pipefail
source "./config.sh"
set +u
eval "$(conda shell.bash hook)"
conda activate amplicon-seq
set -u

# Create each index only if its files are missing or empty
for ext in amb ann bwt pac sa; do
    [[ -s "${REF}.${ext}" ]] || { bwa index "$REF"; break; }
done
[[ -s "${REF}.fai" ]] || samtools faidx "$REF"

echo "Reference indexes ready: $REF"
