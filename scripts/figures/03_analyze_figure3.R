# Analyze the data used by Figure 3 and write its plotting tables.
#
# The Gambia/Kenya/Ghana estimates and counts come from the cited papers.

library(dplyr)

source("scripts/get_pfsa_counts.R")


group_levels <- c(
  "Gambia (severe)", "Kenya (severe)", "Ghana (mild)",
  "Cameroon (asymptomatic)"
)
variant_levels <- c("Pfsa1", "Pfsa3")

# Cameroon estimates used by Fig3.R.
hptest <- read.csv(
  "results/hptest/hbb_vs_Pf_no_covariates.hptest.csv",
  comment.char = "#", check.names = TRUE
)
hptest_required <- c(
  "predictor.position", "outcome.position", "N",
  "gen.beta_1.add.outcome.1", "gen.se_1", "gen.pvalue_1"
)

cameroon_raw <- hptest[
  hptest$predictor.position == 5227002 &
    hptest$outcome.position %in% c(631190, 1058035),
  ,
  drop = FALSE
]
cameroon <- data.frame(
  Variant = ifelse(
    cameroon_raw$outcome.position == 631190, "Pfsa1", "Pfsa3"
  ),
  Group = "Cameroon (asymptomatic)",
  Beta_log_OR = cameroon_raw$gen.beta_1.add.outcome.1,
  Standard_error = cameroon_raw$gen.se_1,
  Odds_ratio = exp(cameroon_raw$gen.beta_1.add.outcome.1),
  CI_95_lower = exp(
    cameroon_raw$gen.beta_1.add.outcome.1 -
      1.96 * cameroon_raw$gen.se_1
  ),
  CI_95_upper = exp(
    cameroon_raw$gen.beta_1.add.outcome.1 +
      1.96 * cameroon_raw$gen.se_1
  ),
  N = cameroon_raw$N,
  P_value = cameroon_raw$gen.pvalue_1,
  Estimate_source = paste(
    "Unadjusted HPTest; OR/CI calculated from beta and standard error",
    "output from HPTest in this study"
  ),
  stringsAsFactors = FALSE
)
cameroon <- cameroon[
  match(variant_levels, cameroon$Variant), , drop = FALSE
]

# Previously published studies used in Figure 3. Gambia and Kenya ORs are
# from Band et al. 2022 Table S2 (additive model), doi:10.1038/s41586-021-04288-3
# Ghana ORs are calculated from the posterior mode (beta) and standard error
# from Hamilton et al. 2026 Supplementary Table 1, doi:10.1186/s12936-026-05956-3

published_plot <- data.frame(
  Variant = rep(c("Pfsa1", "Pfsa3"), each = 3),
  Group = rep(
    c("Gambia (severe)", "Kenya (severe)", "Ghana (mild)"), 2
  ),
  OR = c(3.41, 24.64, NA, 22.49, 18.45, NA),
  CI_lower = c(1.51, 10.66, NA, 3.04, 8.33, NA),
  CI_upper = c(7.72, 59.96, NA, 166.27, 40.89, NA),
  N = c(2025, 1741, 1286, 1958, 1719, 1281),
  P_value = c(2.38e-03, 4.52e-14, 1.46e-11, 1.43e-03, 4.52e-13, 8.14e-16),
  Reported_beta = c(NA, NA, 2.815, NA, NA, 3.167),
  Reported_SE = c(NA, NA, 0.417, NA, NA, 0.393),
  stringsAsFactors = FALSE
)

# Hamilton et al. 2026 reports beta and its standard error for Ghana. Convert
# these log-odds estimates to odds ratios and Wald 95% confidence intervals.
ghana_rows <- published_plot$Group == "Ghana (mild)"
published_plot$OR[ghana_rows] <- round(
  exp(published_plot$Reported_beta[ghana_rows]), 2
)
published_plot$CI_lower[ghana_rows] <- round(
  exp(
    published_plot$Reported_beta[ghana_rows] -
      1.96 * published_plot$Reported_SE[ghana_rows]
  ),
  2
)
published_plot$CI_upper[ghana_rows] <- round(
  exp(
    published_plot$Reported_beta[ghana_rows] +
      1.96 * published_plot$Reported_SE[ghana_rows]
  ),
  2
)

published_plot_rows <- data.frame(
  Variant = published_plot$Variant,
  Group = published_plot$Group,
  Beta_log_OR = ifelse(
    is.na(published_plot$Reported_beta),
    log(published_plot$OR),
    published_plot$Reported_beta
  ),
  Standard_error = published_plot$Reported_SE,
  Odds_ratio = published_plot$OR,
  CI_95_lower = published_plot$CI_lower,
  CI_95_upper = published_plot$CI_upper,
  N = published_plot$N,
  P_value = published_plot$P_value,
  Estimate_source = ifelse(
    published_plot$Group == "Ghana (mild)",
    paste(
      "Calculated OR and CI using published beta and standard error from",
      "Hamilton et al 2026, Supplementary Table 1"
    ),
    paste(
      "Published OR/CI (additive) from Band et al 2022,",
      "Supplementary Table 2"
    )
  ),
  stringsAsFactors = FALSE
)

fig3a <- rbind(published_plot_rows, cameroon)
fig3a <- fig3a[
  order(
    match(fig3a$Variant, variant_levels),
    match(fig3a$Group, group_levels)
  ),
  ,
  drop = FALSE
]
# Gambia and Kenya Pfsa genotype counts by HbS genotype, from Band et al.
# 2022 Table S2.
published_long <- data.frame(
  Country = rep(c("Gambia", "Kenya"), each = 4, times = 2),
  Pfsa = rep(
    c("Pfsa1_plus", "Pfsa1_minus", "Pfsa3_plus", "Pfsa3_minus"), 4
  ),
  HbS_genotype = rep(c("AA", "AS"), each = 8),
  count = c(
    560, 1447, 828, 1111,
    200, 1504, 220, 1463,
    10, 5, 16, 0,
    25, 4, 24, 5
  ),
  Study_source = "Previously published in Band et al 2022",
  stringsAsFactors = FALSE
)

# Cameroon: same counts as Fig2AB, restricted to
# the plus (1/1) / minus (0/0) categories plotted in Figure 3B.
metadata <- read.table(
  "data/metadata.tsv", check.names = FALSE, header = TRUE, sep = "\t"
)

cameroon_counts_all <- rbind(
  get_pfsa_counts(metadata, "Pfsa1", "Pf3D7_02_v3.631190_gt"),
  get_pfsa_counts(metadata, "Pfsa3", "Pf3D7_11_v3.1058035_gt")
)
cameroon_plus_minus <- cameroon_counts_all[
  cameroon_counts_all$Parasite_genotype %in% c("0/0", "1/1"),
  ,
  drop = FALSE
]
cameroon_counts <- data.frame(
  Country = "Cameroon",
  Pfsa = paste0(
    cameroon_plus_minus$Site, "_",
    ifelse(
      cameroon_plus_minus$Parasite_genotype == "1/1", "plus", "minus"
    )
  ),
  HbS_genotype = as.character(cameroon_plus_minus$HbS_genotype),
  count = cameroon_plus_minus$count,
  Study_source = "This study",
  stringsAsFactors = FALSE
)

# Ghana: Pfsa genotype counts are estimated from Hamilton et al. 2026
# Supplementary Figure 3 (Pfsa frequency by sickle status) multiplied by
# the genotyped N per HbS genotype from Supplementary Table 1.

ghana_sample_sizes <- c(
  Pfsa1_HbAA = 1255, Pfsa1_HbAS = 30,
  Pfsa3_HbAA = 1248, Pfsa3_HbAS = 32
)
ghana_pfsa_frequencies <- c(
  Pfsa1_HbAA = 0.03, Pfsa1_HbAS = 0.39,
  Pfsa3_HbAA = 0.03, Pfsa3_HbAS = 0.49
)
ghana_plus_counts <- round(ghana_pfsa_frequencies * ghana_sample_sizes)

ghana_counts <- data.frame(
  Country = "Ghana",
  Pfsa = rep(
    c("Pfsa1_plus", "Pfsa1_minus", "Pfsa3_plus", "Pfsa3_minus"),
    each = 2
  ),
  HbS_genotype = rep(c("AA", "AS"), 4),
  count = unname(c(
    ghana_plus_counts[c("Pfsa1_HbAA", "Pfsa1_HbAS")],
    ghana_sample_sizes[c("Pfsa1_HbAA", "Pfsa1_HbAS")] -
      ghana_plus_counts[c("Pfsa1_HbAA", "Pfsa1_HbAS")],
    ghana_plus_counts[c("Pfsa3_HbAA", "Pfsa3_HbAS")],
    ghana_sample_sizes[c("Pfsa3_HbAA", "Pfsa3_HbAS")] -
      ghana_plus_counts[c("Pfsa3_HbAA", "Pfsa3_HbAS")]
  )),
  Study_source = paste(
    "Estimated from Hamilton et al. 2026 by taking ",
    "Supplementary Figure 3 Pfsa frequency by sickle status,",
    "multiplied by genotyped N from Supplementary Table 1, rounded",
    "to whole individuals so plus + minus equals the genotyped N"
  ),
  stringsAsFactors = FALSE
)

fig3b <- rbind(published_long, ghana_counts, cameroon_counts)
split_pfsa <- strsplit(fig3b$Pfsa, "_", fixed = TRUE)
fig3b$Variant <- vapply(split_pfsa, `[`, character(1), 1)
fig3b$Pfsa_status <- vapply(split_pfsa, `[`, character(1), 2)
fig3b <- fig3b %>%
  group_by(Country, Variant, HbS_genotype) %>%
  mutate(
    total_within_HbS = sum(count),
    proportion = count / total_within_HbS
  ) %>%
  ungroup() %>%
  as.data.frame()
fig3b <- fig3b[, c(
  "Country", "Variant", "HbS_genotype", "Pfsa_status", "count",
  "total_within_HbS", "proportion", "Study_source"
)]

fig3a_output <- fig3a[, c(
  "Variant", "Group", "Odds_ratio", "CI_95_lower", "CI_95_upper",
  "N", "P_value", "Estimate_source"
)]
dir.create("results/figures/data", recursive = TRUE, showWarnings = FALSE)
write.table(
  fig3a_output, "results/figures/data/Figure3A_odds_ratios.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  fig3b, "results/figures/data/Figure3B_genotype_counts.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)
