# Pfsa and sickle hemoglobin in asymptomatic *P.falciparum* infection

Scripts for figures and analysis in the manuscript:

Hopson HD, et al. *Sickle cell status skews malaria parasite genotype at infection*. bioRxiv (2025). doi: https://doi.org/10.1101/2025.09.09.675015

## Requirements
- R (tested with 4.4.1) and the following packages: `ggplot2`, `dplyr`, `patchwork`, `cowplot`, `scales`, `tidyr`

## Organization

```
scripts/       one script per figure: builds TSV data tables, then plots them
figure_data/   one tab-separated data file per figure panel or panel group
figures/       generated panel PDFs 
inputs/        the external files the scripts read to build the data tables

```

## Inputs

- `inputs/metadata.tsv` - reduced export of the per-sample metadata from Supplementary Table 9 containing only columns used by figure scripts

- `inputs/hptest_no_covariates.csv` — unadjusted association-test output, used for
  the Cameroon Pfsa1/Pfsa3 odds ratios in Figure 3.
