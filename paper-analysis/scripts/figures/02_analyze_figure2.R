# Analyze the data used by Figure 2 and write its plotting tables.

library(dplyr)

source("scripts/get_pfsa_counts.R")


metadata <- read.table(
  "data/analysis_metadata.tsv", check.names = FALSE, header = TRUE, sep = "\t"
)

fig2ab <- rbind(
  get_pfsa_counts(metadata, "Pfsa1", "Pf3D7_02_v3.631190_gt"),
  get_pfsa_counts(metadata, "Pfsa3", "Pf3D7_11_v3.1058035_gt")
)

pfsa1_mixed <- metadata$HbS_gt %in% c("AA", "AS") &
  metadata$Pf3D7_02_v3.631190_gt == "0/1" &
  !is.na(metadata$pfsa1_vaf)
pfsa3_mixed <- metadata$HbS_gt %in% c("AA", "AS") &
  metadata$Pf3D7_11_v3.1058035_gt == "0/1" &
  !is.na(metadata$pfsa3_vaf)

fig2cd <- rbind(
  data.frame(
    Site = "Pfsa1",
    Sequencing_sample_ID = metadata$seq_sample_id[pfsa1_mixed],
    HbS_genotype = metadata$HbS_gt[pfsa1_mixed],
    VAF = metadata$pfsa1_vaf[pfsa1_mixed]
  ),
  data.frame(
    Site = "Pfsa3",
    Sequencing_sample_ID = metadata$seq_sample_id[pfsa3_mixed],
    HbS_genotype = metadata$HbS_gt[pfsa3_mixed],
    VAF = metadata$pfsa3_vaf[pfsa3_mixed]
  )
)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.table(
  fig2ab, "results/tables/Figure2AB_genotype_counts.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  fig2cd, "results/tables/Figure2CD_VAF.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)
