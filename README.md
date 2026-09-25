# Amplicon sequencing pipeline and host–parasite association analysis

Code accompanying: Hopson HD, et al. “Sickle cell status skews malaria parasite genotype at
infection.” *bioRxiv* (2025). [doi:10.1101/2025.09.09.675015](https://doi.org/10.1101/2025.09.09.675015)

Amplicon sequencing of 2,246 dried blood spots from asymptomatic children in Mfou, Cameroon, 
targeting human variants associated with malaria (including the sickle cell allele) and *P. falciparum* loci including the sickle-associated (Pfsa) alleles. 

## Overview

```text
amplicon_pipeline/    FASTQ processing, QC, infection calling, variant calling (human and P. falciparum)
analysis/             Manuscript analyses and figures
```

See [`amplicon_pipeline/README.md`](amplicon_pipeline/README.md) and [`analysis/README.md`](analysis/README.md)

The variant calls output from `amplicon_pipeline/` are used in `analysis/`: 

- Human: `amplicon_pipeline/output/human/<amplicon>.human.filtered.final.vcf.gz`
- *P. falciparum*: `amplicon_pipeline/output/parasite/pf.all.filtered.final.vcf.gz`

## Data availability

Sequencing reads are deposited in NCBI SRA under BioProject
[PRJNA1503656](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1503656)