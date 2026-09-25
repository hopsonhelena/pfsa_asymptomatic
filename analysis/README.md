# Manuscript analysis

Code and data used for the manuscript figures and association analyses.

## Overview

```text
data/                 analysis inputs (see its README)
scripts/
  figures/            analysis and plotting scripts
  association/        association test scripts
  *.R                 analysis functions
results/
  figures/            generated PDFs
  tables/             generated source-data tables
  association/
    hptest/           genotype x genotype association outputs
    snptest/          genotype x phenotype association outputs
```

## Requirements

Figures were created with R 4.4.1 and require `ggplot2`, `dplyr`, `patchwork`,
`cowplot`, `scales`, and `tidyr`.

The association analyses require [SNPTEST 2.5.2](https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest) and HPTEST 2.2.0, distributed as part of [QCTOOL](https://github.com/gavinband/qctool).

## Association tests

From this directory (`analysis/`), run:

```bash
bash scripts/association/run_hptest.sh
bash scripts/association/run_snptest_1701_infected_samples.sh
bash scripts/association/run_snptest_2246_samples.sh
```


## Reproduce the figures

Run each figure's analysis script followed by its plotting script from this
directory. For example:

```bash
Rscript scripts/figures/01_analyze_figure1.R
Rscript scripts/figures/01_plot_figure1.R
```

Generated tables are written to `results/tables/` and plot PDFs to
`results/figures/`.

Figure 4 also writes `results/tables/Figure4A_LD_r2.tsv`, containing
Pfsa1/Pfsa3 linkage disequilibrium statistics.

Input files are documented in [`data/README.md`](data/README.md).
