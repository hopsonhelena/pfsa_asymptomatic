# Sickle cell status and malaria parasite genotype in asymptomatic infections

Code for processing and analyzing amplicon sequencing data from Hopson, H.D., Herbert-Mainero, A.,
Bouopda-Tuedom, A.G. et al. “Sickle cell haemoglobin status shapes malaria parasite genotype in
asymptomatic infections.” *Nat Microbiol* (2026).
[doi:10.1038/s41564-026-02502-4](https://doi.org/10.1038/s41564-026-02502-4)

Amplicon sequencing of 2,246 dried blood spots from asymptomatic children in Mfou, Cameroon, 
targeting human variants associated with malaria (including the sickle cell allele) and *P. falciparum* loci including the sickle-associated (Pfsa) alleles. 

## Overview

```text
amplicon-seq/         FASTQ processing, QC, infection calling, variant calling (human and P. falciparum)
paper-analysis/       Manuscript analyses and figures
```

See [`amplicon-seq/README.md`](amplicon-seq/README.md) and [`paper-analysis/README.md`](paper-analysis/README.md)

The variant calls output from `amplicon-seq/` are used in `paper-analysis/`: 

- Human: `amplicon-seq/output/human/<amplicon>.human.filtered.final.vcf.gz`
- *P. falciparum*: `amplicon-seq/output/parasite/pf.all.filtered.final.vcf.gz`

## Data availability

Sequencing reads are deposited in NCBI SRA under BioProject
[PRJNA1503656](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1503656)