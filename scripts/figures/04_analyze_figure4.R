# Analyze the data used by Figure 4 and write its plotting and LD tables.

source("scripts/calculate_observed_expected_ld.R")
source("scripts/calculate_ld_r2.R")

metadata <- read.table(
  "data/metadata.tsv", check.names = FALSE, header = TRUE, sep = "\t"
)
pfsa1_column <- "Pf3D7_02_v3.631190_gt"
pfsa3_column <- "Pf3D7_11_v3.1058035_gt"

observed_expected <- rbind(
  calculate_observed_expected_ld(
    data = metadata,
    all_sample_data = metadata,
    group = "All",
    pfsa1_column = pfsa1_column,
    pfsa3_column = pfsa3_column
  ),
  calculate_observed_expected_ld(
    data = metadata[metadata$HbS_gt == "AA", ],
    all_sample_data = metadata,
    group = "AA",
    pfsa1_column = pfsa1_column,
    pfsa3_column = pfsa3_column
  ),
  calculate_observed_expected_ld(
    data = metadata[metadata$HbS_gt == "AS", ],
    all_sample_data = metadata,
    group = "AS",
    pfsa1_column = pfsa1_column,
    pfsa3_column = pfsa3_column
  )
)

ld_statistics <- rbind(
  calculate_ld_r2(metadata, "All", pfsa1_column, pfsa3_column),
  calculate_ld_r2(
    metadata[metadata$HbS_gt == "AA", ],
    "HbAA", pfsa1_column, pfsa3_column
  ),
  calculate_ld_r2(
    metadata[metadata$HbS_gt == "AS", ],
    "HbAS", pfsa1_column, pfsa3_column
  )
)

mixed_either <- metadata[
  metadata$Infected == 1 &
    (
      metadata[[pfsa1_column]] == "0/1" |
        metadata[[pfsa3_column]] == "0/1"
    ) &
    metadata$HbS_gt %in% c("AA", "AS") &
    complete.cases(metadata[, c("pfsa1_vaf", "pfsa3_vaf")]),
  ,
  drop = FALSE
]
paired_vaf <- mixed_either[, c(
  "Seq.Sample.ID", "Sample.ID", "HbS_gt", "pfsa1_vaf", "pfsa3_vaf"
)]
names(paired_vaf) <- c(
  "Sequencing_sample_ID", "Study_sample_ID", "HbS_genotype",
  "Pfsa1_VAF", "Pfsa3_VAF"
)

dir.create("results/figures/data", recursive = TRUE, showWarnings = FALSE)
write.table(
  observed_expected, "results/figures/data/Figure4A_observed_expected.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  ld_statistics, "results/figures/data/Figure4A_LD_r2.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  paired_vaf, "results/figures/data/Figure4B_paired_VAF.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)
