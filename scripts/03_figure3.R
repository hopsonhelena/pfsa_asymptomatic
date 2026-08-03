# Create and plot data for main Figure 3.
#
# The Gambia/Kenya/Ghana odds ratios and genotype counts are from the published papers 
# cited in the relevant comments.

library(ggplot2)
library(patchwork)
library(dplyr)


# ---- Create table for Figure 3 ----

group_levels <- c(
  "Gambia (severe)", "Kenya (severe)", "Ghana (mild)",
  "Cameroon (asymptomatic)"
)
variant_levels <- c("Pfsa1", "Pfsa3")

# Cameroon estimates used by Fig3.R.
hptest <- read.csv(
  "inputs/hptest_no_covariates.csv",
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
  "inputs/metadata.tsv", check.names = FALSE, header = TRUE, sep = "\t"
)

get_pfsa_counts <- function(metadata, site, genotype_column) {
  keep <- metadata$HbS_gt %in% c("AA", "AS") &
    !is.na(metadata[[genotype_column]]) &
    metadata[[genotype_column]] != "./."
  selected <- metadata[
    keep, c("HbS_gt", genotype_column), drop = FALSE
  ]
  names(selected) <- c("HbS_genotype", "Parasite_genotype")
  result <- as.data.frame(
    table(
      HbS_genotype = factor(selected$HbS_genotype, c("AA", "AS")),
      Parasite_genotype = factor(
        selected$Parasite_genotype, c("0/0", "0/1", "1/1")
      )
    ),
    responseName = "count"
  )
  result <- result[result$count > 0, , drop = FALSE]
  result <- result %>%
    group_by(HbS_genotype) %>%
    mutate(
      total_within_HbS = sum(count),
      proportion = count / total_within_HbS
    ) %>%
    ungroup() %>%
    as.data.frame()
  result$Site <- site
  result[, c(
    "Site", "HbS_genotype", "Parasite_genotype",
    "count", "total_within_HbS", "proportion"
  )]
}

get_cameroon_counts <- function(site, genotype_column) {
  fig2_counts <- get_pfsa_counts(metadata, site, genotype_column)
  plus_minus <- fig2_counts[
    fig2_counts$Parasite_genotype %in% c("0/0", "1/1"), , drop = FALSE
  ]
  data.frame(
    Country = "Cameroon",
    Pfsa = paste0(
      site, "_",
      ifelse(plus_minus$Parasite_genotype == "1/1", "plus", "minus")
    ),
    HbS_genotype = as.character(plus_minus$HbS_genotype),
    count = plus_minus$count,
    Study_source = "This study",
    stringsAsFactors = FALSE
  )
}
cameroon_counts <- rbind(
  get_cameroon_counts("Pfsa1", "Pf3D7_02_v3.631190_gt"),
  get_cameroon_counts("Pfsa3", "Pf3D7_11_v3.1058035_gt")
)

# Ghana: Pfsa genotype counts are estimated from Hamilton et al. 2026 
# Supplementary Figure 3 (Pfsa frequency by sickle status) multiplied by 
# the genotyped N per HbS genotype from Supplementary Table 1.

get_ghana_counts <- function(
  site, sickle_freq, nonsickle_freq, HbAA_n, HbAS_n
) {
  HbAA_plus <- round(nonsickle_freq * HbAA_n)
  HbAS_plus <- round(sickle_freq * HbAS_n)
  data.frame(
    Country = "Ghana",
    Pfsa = paste0(site, "_", c("plus", "plus", "minus", "minus")),
    HbS_genotype = c("AA", "AS", "AA", "AS"),
    count = c(
      HbAA_plus, HbAS_plus,
      HbAA_n - HbAA_plus, HbAS_n - HbAS_plus
    ),
    Study_source = paste(
      "Estimated from Hamilton et al. 2026 by taking ",
      "Supplementary Figure 3 Pfsa frequency by sickle status,",
      "multiplied by genotyped N from Supplementary Table 1, rounded",
      "to whole individuals so plus + minus equals the genotyped N"
    ),
    stringsAsFactors = FALSE
  )
}
ghana_counts <- rbind(
  get_ghana_counts(
    "Pfsa1",
    sickle_freq = 0.39, nonsickle_freq = 0.03,
    HbAA_n = 1255, HbAS_n = 30
  ),
  get_ghana_counts(
    "Pfsa3",
    sickle_freq = 0.49, nonsickle_freq = 0.03,
    HbAA_n = 1248, HbAS_n = 32
  )
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
write.table(
  fig3a_output, "figure_data/Figure3A_odds_ratios.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  fig3b, "figure_data/Figure3B_genotype_counts.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

# ---- Plot figure ----

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
theme_set(theme_minimal(base_size = 7))

or_data <- fig3a_output
count_data <- fig3b
country_levels <- c("Gambia", "Kenya", "Ghana", "Cameroon")
or_data$Variant <- factor(or_data$Variant, variant_levels)
or_data$Group <- factor(or_data$Group, group_levels)
count_data$Country <- factor(count_data$Country, country_levels)
count_data$HbS_genotype <- factor(
  count_data$HbS_genotype, c("AA", "AS")
)

make_or_plot <- function(data, variant_name) {
  plot_data <- data[data$Variant == variant_name, , drop = FALSE]

  ggplot(
    plot_data,
    aes(
      y = Group, x = Odds_ratio,
      xmin = CI_95_lower, xmax = CI_95_upper, color = Group
    )
  ) +
    geom_pointrange(linewidth = 0.5, size = 0.3) +
    annotate(
      "segment",
      x = 0.5, xend = 600,
      y = 3.5, yend = 3.5,
      color = "grey60",
      linewidth = 0.2,
      linetype = "dashed"
    ) +
    geom_vline(xintercept = 1, color = "black", linewidth = 0.2) +
    scale_x_log10(
      breaks = c(0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512),
      limits = c(0.5, 600),
      oob = scales::oob_keep
    ) +
    scale_color_manual(values = c(
      "Cameroon (asymptomatic)" = "black",
      "Ghana (mild)" = "grey60",
      "Kenya (severe)" = "grey60",
      "Gambia (severe)" = "grey60"
    )) +
    scale_y_discrete(labels = rep("", length(group_levels))) +
    coord_cartesian(clip = "off") +
    labs(x = "OR", y = NULL, title = variant_name) +
    theme_minimal(base_size = 7) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 7, hjust = 0.5),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      plot.margin = margin(5, 5, 5, -10)
    )
}

pfsa1_or_plot <- make_or_plot(or_data, "Pfsa1")
pfsa3_or_plot <- make_or_plot(or_data, "Pfsa3")

count_labels <- unique(
  count_data[, c(
    "Country", "Variant", "HbS_genotype", "total_within_HbS"
  )]
)
make_count_plot <- function(
  data, labels, variant_name, show_legend = FALSE
) {
  plot_data <- data[data$Variant == variant_name, , drop = FALSE]
  plot_labels <- labels[labels$Variant == variant_name, , drop = FALSE]

  plot <- ggplot(
    plot_data,
    aes(x = proportion, y = Country, fill = Pfsa_status)
  ) +
    geom_col(color = "black", linewidth = 0.2, width = 0.7) +
    annotate(
      "segment",
      x = 0, xend = 1.1,
      y = 3.5, yend = 3.5,
      color = "grey60",
      linewidth = 0.2,
      linetype = "dashed"
    ) +
    geom_text(
      data = plot_labels,
      aes(
        x = 1.02, y = Country, label = total_within_HbS
      ),
      inherit.aes = FALSE,
      hjust = 0,
      size = 2.2
    ) +
    facet_grid(~HbS_genotype) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(
      breaks = seq(0, 1, 0.25),
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0.15, 0.1))
    ) +
    scale_y_discrete(labels = rep("", length(country_levels))) +
    scale_fill_manual(
      name = NULL,
      values = c("minus" = "#009E73", "plus" = "#D55E00"),
      labels = c("minus" = "Pfsa-", "plus" = "Pfsa+")
    ) +
    labs(x = "Proportion", y = NULL, fill = NULL) +
    theme_minimal(base_size = 7) +
    theme(
      panel.spacing.x = grid::unit(0, "lines"),
      strip.text.x = element_text(hjust = 0.5),
      strip.text = element_text(size = 7),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(10, 30, 5, 5)
    )

  if (show_legend) {
    plot + theme(
      legend.key.size = grid::unit(3, "mm"),
      legend.key.width = grid::unit(3, "mm"),
      legend.box.spacing = grid::unit(3, "mm")
    )
  } else {
    plot + theme(legend.position = "none")
  }
}

pfsa1_count_plot <- make_count_plot(
  count_data, count_labels, "Pfsa1", show_legend = FALSE
)
pfsa3_count_plot <- make_count_plot(
  count_data, count_labels, "Pfsa3", show_legend = TRUE
)

study_labels <- ggplot() +
  annotate(
    "text", x = 0.28, y = 4.3,
    label = "Cameroon", hjust = 0,
    size = 2.2, color = "black", fontface = "bold"
  ) +
  annotate(
    "text", x = 0.28, y = 4.1,
    label = "(asymptomatic)", hjust = 0,
    size = 2.2, color = "black", fontface = "bold"
  ) +
  annotate(
    "text", x = 0.28, y = 3.5,
    label = "previous reports", hjust = 0,
    size = 2.2, color = "grey30", fontface = "italic"
  ) +
  annotate(
    "text", x = 0.31, y = 2.8,
    label = "Ghana (mild)", hjust = 0,
    size = 2.2, color = "grey30"
  ) +
  annotate(
    "text", x = 0.31, y = 1.8,
    label = "Kenya (severe)", hjust = 0,
    size = 2.2, color = "grey30"
  ) +
  annotate(
    "text", x = 0.31, y = 0.8,
    label = "Gambia (severe)", hjust = 0,
    size = 2.2, color = "grey30"
  ) +
  coord_cartesian(
    xlim = c(0, 0.5),
    ylim = c(0.3, 4.5),
    clip = "off"
  ) +
  theme_void() +
  theme(plot.margin = margin(5, 0, 5, 5))

top_row <- (study_labels | pfsa1_or_plot | pfsa3_or_plot) +
  plot_layout(widths = c(1.1, 1.8, 1.8))
bottom_row <- (
  study_labels | pfsa1_count_plot | pfsa3_count_plot
) +
  plot_layout(widths = c(1.1, 1.8, 1.8))

ggsave(
  file.path("figures", "Figure3A.pdf"),
  plot = top_row, device = "pdf",
  width = 178, height = 50, units = "mm", bg = "white"
)
ggsave(
  file.path("figures", "Figure3B.pdf"),
  plot = bottom_row, device = "pdf",
  width = 178, height = 55, units = "mm", bg = "white"
)
