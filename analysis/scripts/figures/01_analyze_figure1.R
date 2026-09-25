# Analyze the data used by Figure 1 and write its plotting tables.

library(dplyr)

source("scripts/count_rate_rows.R")


metadata <- read.table(
  "data/analysis_metadata.tsv", check.names = FALSE, header = TRUE, sep = "\t"
)
metadata_aa_as <- metadata[
  metadata$HbS_gt %in% c("AA", "AS"), , drop = FALSE
]

infection <- count_rate_rows(metadata_aa_as, "HbS_gt", "Infected")
gametocytes <- count_rate_rows(
  metadata_aa_as[metadata_aa_as$Infected == 1, , drop = FALSE],
  "HbS_gt", "Gametocytes_YN"
)
fig1a <- merge(
  infection, gametocytes, by = "HbS_genotype",
  suffixes = c("_infection", "_gametocyte"), sort = FALSE
)
fig1a <- fig1a[match(c("AA", "AS"), fig1a$HbS_genotype), ]
fig1a <- data.frame(
  HbS_genotype = fig1a$HbS_genotype,
  N = fig1a$N_infection,
  N_infected = fig1a$count_infection,
  proportion_infected = fig1a$proportion_infection,
  CI_95_lower = fig1a$CI_95_lower_infection,
  CI_95_upper = fig1a$CI_95_upper_infection,
  N_infected_gametocytes = fig1a$count_gametocyte,
  proportion_infected_gametocytes = fig1a$proportion_gametocyte,
  check.names = FALSE
)

ama1_evaluable <- metadata[
  metadata$HbS_gt %in% c("AA", "AS") & !is.na(metadata$COI_ama1),
  ,
  drop = FALSE
]
ama1 <- count_rate_rows(
  ama1_evaluable, "HbS_gt", "monopolyCOI_ama1"
)
names(ama1)[names(ama1) == "count"] <- "polyclonal_count"
names(ama1)[names(ama1) == "count_absent"] <- "monoclonal_count"
names(ama1)[names(ama1) == "proportion"] <- "polyclonal_proportion"
ama1$Amplicon <- "AMA1"

sera2_evaluable <- metadata[
  metadata$HbS_gt %in% c("AA", "AS") & !is.na(metadata$COI_sera2),
  ,
  drop = FALSE
]
sera2 <- count_rate_rows(
  sera2_evaluable, "HbS_gt", "monopolyCOI_sera2"
)
names(sera2)[names(sera2) == "count"] <- "polyclonal_count"
names(sera2)[names(sera2) == "count_absent"] <- "monoclonal_count"
names(sera2)[names(sera2) == "proportion"] <- "polyclonal_proportion"
sera2$Amplicon <- "SERA2"

fig1b_columns <- c(
  "Amplicon", "HbS_genotype", "N", "polyclonal_count",
  "monoclonal_count", "polyclonal_proportion",
  "CI_95_lower", "CI_95_upper"
)
fig1b <- rbind(ama1[fig1b_columns], sera2[fig1b_columns])

metadata$Parasitemia_plot <- ifelse(
  metadata$Parasitemia %in% c("Pf+++++", "Pf++++++"),
  "Pf++++", metadata$Parasitemia
)
fig1c_metadata <- metadata[
  metadata$HbS_gt %in% c("AA", "AS"), , drop = FALSE
]
fig1c <- as.data.frame(
  table(
    HbS_genotype = factor(fig1c_metadata$HbS_gt, c("AA", "AS")),
    Parasitemia_category = factor(
      fig1c_metadata$Parasitemia_plot,
      c("0", "Pf+", "Pf++", "Pf+++", "Pf++++")
    )
  ),
  responseName = "count"
)
fig1c <- fig1c[fig1c$count > 0, , drop = FALSE]
fig1c <- fig1c %>%
  group_by(HbS_genotype) %>%
  mutate(
    total_within_HbS = sum(count),
    proportion = count / total_within_HbS
  ) %>%
  ungroup() %>%
  as.data.frame()

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.table(
  fig1a, "results/tables/Figure1A_rates.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  fig1b, "results/tables/Figure1B_polyclonality.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  fig1c, "results/tables/Figure1C_parasitemia.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)
