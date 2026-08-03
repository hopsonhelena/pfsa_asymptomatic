# Create and plot data for main Figure 1.

library(ggplot2)
library(dplyr)


# ---- Create table for Figure 1 ----

metadata <- read.table(
  "inputs/metadata.tsv", check.names = FALSE, header = TRUE, sep = "\t"
)
metadata_aa_as <- metadata[
  metadata$HbS_gt %in% c("AA", "AS"), , drop = FALSE
]

count_rate_rows <- function(data, group_column, outcome_column) {
  groups <- split(data, data[[group_column]], drop = TRUE)
  output <- lapply(names(groups), function(group) {
    values <- groups[[group]][[outcome_column]]
    n <- length(values)
    successes <- sum(values == 1, na.rm = TRUE)
    estimate <- successes / n
    se <- sqrt(estimate * (1 - estimate) / n)
    data.frame(
      HbS_genotype = group,
      N = n,
      count = successes,
      count_absent = n - successes,
      proportion = estimate,
      CI_95_lower = max(0, estimate - 1.96 * se),
      CI_95_upper = min(1, estimate + 1.96 * se)
    )
  })
  do.call(rbind, output)
}

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

get_polyclonal <- function(amplicon, coi_column, outcome_column) {
  evaluable <- metadata[
    metadata$HbS_gt %in% c("AA", "AS") &
      !is.na(metadata[[coi_column]]),
    ,
    drop = FALSE
  ]
  result <- count_rate_rows(evaluable, "HbS_gt", outcome_column)
  names(result)[names(result) == "count"] <- "polyclonal_count"
  names(result)[names(result) == "count_absent"] <- "monoclonal_count"
  names(result)[names(result) == "proportion"] <- "polyclonal_proportion"
  result$Amplicon <- amplicon
  result[, c(
    "Amplicon", "HbS_genotype", "N", "polyclonal_count",
    "monoclonal_count", "polyclonal_proportion",
    "CI_95_lower", "CI_95_upper"
  )]
}
fig1b <- rbind(
  get_polyclonal("AMA1", "COI_ama1", "monopolyCOI_ama1"),
  get_polyclonal("SERA2", "COI_sera2", "monopolyCOI_sera2")
)

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

write.table(
  fig1a, "figure_data/Figure1A_rates.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  fig1b, "figure_data/Figure1B_polyclonality.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  fig1c, "figure_data/Figure1C_parasitemia.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

# ---- Plot figure ----

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
theme_set(theme_minimal(base_size = 7))

fig1a$HbS_genotype <- factor(fig1a$HbS_genotype, c("AA", "AS"))
stacked_df <- fig1a %>%
  transmute(
    HbS_genotype,
    gametocyte = proportion_infected_gametocytes,
    remainder = pmax(
      proportion_infected - proportion_infected_gametocytes, 0
    )
  ) %>%
  tidyr::pivot_longer(
    cols = c(gametocyte, remainder),
    names_to = "segment",
    values_to = "rate"
  )
stacked_df$segment <- factor(
  stacked_df$segment, levels = c("gametocyte", "remainder")
)

aplot <- ggplot() +
  geom_col(
    data = stacked_df,
    aes(x = HbS_genotype, y = rate, fill = segment),
    width = 0.6,
    alpha = 0.4,
    color = "black",
    linewidth = 0.2
  ) +
  geom_line(
    data = fig1a,
    aes(
      x = HbS_genotype, y = proportion_infected, group = 1
    ),
    linetype = "dashed"
  ) +
  geom_point(
    data = fig1a,
    aes(x = HbS_genotype, y = proportion_infected),
    shape = 21, fill = "black", size = 1
  ) +
  geom_errorbar(
    data = fig1a,
    aes(
      x = HbS_genotype, ymin = CI_95_lower, ymax = CI_95_upper
    ),
    width = 0.1
  ) +
  scale_fill_manual(
    values = c("gametocyte" = "forestgreen", "remainder" = "gray80"),
    breaks = c("gametocyte", "remainder"),
    labels = c("Present", "Absent"),
    name = "Gametocytes"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(y = "Infection Rate", x = NULL) +
  theme_minimal(base_size = 7) +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.position = "right",
    legend.key.height = grid::unit(0.2, "cm"),
    legend.key.width = grid::unit(0.2, "cm"),
    legend.spacing.x = grid::unit(0.05, "cm"),
    legend.spacing.y = grid::unit(0.05, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 7),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3)))

fig1b$HbS_genotype <- factor(fig1b$HbS_genotype, c("AA", "AS"))
fig1b$Amplicon <- factor(fig1b$Amplicon, c("AMA1", "SERA2"))
bplot <- ggplot(
  fig1b,
  aes(x = HbS_genotype, y = polyclonal_proportion)
) +
  geom_point(shape = 1, position = position_dodge(width = 0.2), size = 1) +
  geom_line(aes(group = 1)) +
  geom_errorbar(
    aes(ymin = CI_95_lower, ymax = CI_95_upper),
    width = 0.2,
    position = position_dodge(width = 0.2)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL, y = "Polyclonal Rate") +
  theme_minimal(base_size = 7) +
  facet_wrap(~Amplicon) +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    strip.text = element_text(size = 7)
  )

fig1c$HbS_genotype <- factor(fig1c$HbS_genotype, c("AA", "AS"))
fig1c$Parasitemia_category <- factor(
  fig1c$Parasitemia_category,
  c("0", "Pf+", "Pf++", "Pf+++", "Pf++++")
)
cplot <- ggplot(
  fig1c,
  aes(x = HbS_genotype, y = proportion, fill = Parasitemia_category)
) +
  geom_col(
    position = position_stack(reverse = TRUE),
    color = "black",
    width = 0.6,
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = c(
      "0" = "white",
      "Pf+" = "#A6BDDB",
      "Pf++" = "#3690C0",
      "Pf+++" = "#0570B0",
      "Pf++++" = "#023858"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL, y = "Proportion", fill = "Parasitemia") +
  theme_minimal(base_size = 7) +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.position = "right",
    legend.key.height = grid::unit(0.2, "cm"),
    legend.key.width = grid::unit(0.2, "cm"),
    legend.spacing.x = grid::unit(0.05, "cm"),
    legend.spacing.y = grid::unit(0.05, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 7),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3)))

ggsave(
  file.path("figures", "Figure1A.pdf"),
  plot = aplot, device = "pdf",
  width = 60, height = 55, units = "mm", bg = "white"
)
ggsave(
  file.path("figures", "Figure1B.pdf"),
  plot = bplot, device = "pdf",
  width = 60, height = 50, units = "mm", bg = "white"
)
ggsave(
  file.path("figures", "Figure1C.pdf"),
  plot = cplot, device = "pdf",
  width = 60, height = 55, units = "mm", bg = "white"
)
