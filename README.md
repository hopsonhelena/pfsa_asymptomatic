# Pfsa and sickle hemoglobin in asymptomatic *P.falciparum* infection

Scripts for the four main figures in the manuscript.

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

| Figure | Data tables | Script | PDF panels |
| --- | --- | --- | --- |
| 1 | `Figure1A_rates.tsv`, `Figure1B_polyclonality.tsv`, `Figure1C_parasitemia.tsv` | `scripts/01_figure1.R` | A-C |
| 2 | `Figure2AB_genotype_counts.tsv`, `Figure2CD_VAF.tsv` | `scripts/02_figure2.R` | A-D |
| 3 | `Figure3A_odds_ratios.tsv`, `Figure3B_genotype_counts.tsv` | `scripts/03_figure3.R` | A-B |
| 4 | `Figure4A_observed_expected.tsv`, `Figure4B_paired_VAF.tsv` | `scripts/04_figure4.R` | A-B |

## Inputs

- `inputs/metadata.tsv` - reduced export of the per-sample metadata from Supplementary Table 9 containing only columns used by figure scripts

- `hptest_no_covariates.csv` — unadjusted association-test output, used for
  the Cameroon Pfsa1/Pfsa3 odds ratios in Figure 3.