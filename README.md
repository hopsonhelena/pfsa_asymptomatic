# Amplicon sequencing pipeline and host–parasite association analysis

This repository contains the analysis workflow for the amplicon sequencing study described in:

Hopson HD, et al. “Sickle cell status skews malaria parasite genotype at
infection.” *bioRxiv* (2025). [doi:10.1101/2025.09.09.675015](https://doi.org/10.1101/2025.09.09.675015)

It is split into two parts: (1) amplicon sequencing pipeline (2) analysis and main figures:

```text
amplicon_pipeline/    FASTQ processing, QC, infection calling, and human and *Plasmodium falciparum* variant calling
analysis/             Manuscript analyses and figures
```

## Run the sequencing pipeline

The pipeline is in `amplicon_pipeline/`. Its configuration,
dependencies, input format, Slurm commands, and output files are documented in
[`amplicon_pipeline/README.md`](amplicon_pipeline/README.md). Briefly:

```bash
cd amplicon_pipeline
conda env create -f environment.yml
# Edit config.sh
source config.sh
bash submit_alignment.sh
```

The final pipeline products are variant calls in human and *P. falciparum*
`amplicon_pipeline/output/human/<amplicon>.human.filtered.final.vcf.gz` and
`amplicon_pipeline/output/parasite/pf.all.filtered.final.vcf.gz`.

## Run the analysis

Requirements, input data, and figure-reproduction instructions are documented
in [`analysis/README.md`](analysis/README.md).
