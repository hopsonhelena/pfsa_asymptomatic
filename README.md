# Pfsa and sickle hemoglobin in asymptomatic *P. falciparum* infection

Code and data supporting the figures and association analyses in:

Hopson HD, et al. “Sickle cell status skews malaria parasite genotype at
infection.” *bioRxiv* (2025). [doi:10.1101/2025.09.09.675015](https://doi.org/10.1101/2025.09.09.675015)

## Repository organization

```text
data/          analysis inputs: metadata, genotypes, and sample files
scripts/
  figures/       paired analysis and plotting scripts for Figures 1-4
  *.R            reusable R analysis functions
  association/   HPTEST and SNPTEST association scripts
results/       generated figures and association outputs
  figures/data/   generated tab-separated source-data tables
  figures/plots/  generated PDF figures and panels
```


## Requirements

Figures were tested with R 4.4.1 and require `ggplot2`, `dplyr`, `patchwork`,
`cowplot`, `scales`, and `tidyr`.

The association analyses require:

- HPTEST 2.2.0, distributed with [QCTOOL](https://www.well.ox.ac.uk/~gav/qctool_v2/)
- SNPTEST 2.5.2


## Reproduce the figures

From the repository root, run each figure's analysis script followed by its
plotting script. For example, to rebuild Figure 1:

```bash
Rscript scripts/figures/01_analyze_figure1.R
Rscript scripts/figures/01_plot_figure1.R
```

Repeat with the corresponding `02`, `03`, and `04` script pairs to rebuild the
other figures. Analysis tables are written to `results/figures/data/`, and PDFs
are written to `results/figures/plots/`.

Figure 4 also writes `results/figures/data/Figure4A_LD_r2.tsv`, containing
Pfsa1/Pfsa3 linkage disequilibrium statistics.


## Input data

- `data/metadata.tsv` is a reduced Supplementary Table 9 containing
  only columns used by the figure scripts.
- The VCF, `.sample`, and text files in `data/` are inputs to the HPTEST and
  SNPTEST scripts.

Figure 3 reads the Cameroon Pfsa1/Pfsa3 estimates directly from
`results/hptest/hbb_vs_Pf_no_covariates.hptest.csv`.
