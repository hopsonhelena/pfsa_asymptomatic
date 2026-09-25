# Amplicon sequencing pipeline and host–parasite association analysis

Code accompanying the amplicon sequencing study described in:

Hopson HD, et al. “Sickle cell status skews malaria parasite genotype at
infection.” *bioRxiv* (2025). [doi:10.1101/2025.09.09.675015](https://doi.org/10.1101/2025.09.09.675015)

## Overview

```text
amplicon_pipeline/   FASTQ processing, QC, infection calling, and human and P. falciparum variant calling
analysis/             Manuscript analyses and figures
```

## Sequencing pipeline

Configuration, dependencies, inputs, Slurm commands, and outputs are
documented in [`amplicon_pipeline/README.md`](amplicon_pipeline/README.md).
Briefly:

```bash
cd amplicon_pipeline
conda env create -f environment.yml
# Edit config.sh for your paths and cluster settings
source ./config.sh
bash submit_alignment.sh
```

Final output variant calls: 

- Human: `amplicon_pipeline/output/human/<amplicon>.human.filtered.final.vcf.gz`
- *P. falciparum*: `amplicon_pipeline/output/parasite/pf.all.filtered.final.vcf.gz`

## Analysis and figures

Requirements, input data, and instructions for reproducing figures are documented
in [`analysis/README.md`](analysis/README.md).
