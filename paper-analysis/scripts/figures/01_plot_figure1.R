# Plot main Figure 1 from its analysis tables.

library(ggplot2)
library(dplyr)
library(patchwork)

fig1a <- read.delim("results/tables/Figure1A_rates.tsv")
fig1b <- read.delim("results/tables/Figure1B_polyclonality.tsv")
fig1c <- read.delim("results/tables/Figure1C_parasitemia.tsv")

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
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, 1)) +
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

figure1 <- (aplot + bplot + cplot) +
  plot_layout(axis_titles = "collect", widths = c(1, 1.3, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 7))

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave(
  file.path("results/figures", "Figure1.pdf"),
  plot = figure1, device = "pdf",
  width = 178, height = 65, units = "mm", bg = "white"
)
ggsave(
  file.path("results/figures", "Figure1A.pdf"),
  plot = aplot, device = "pdf",
  width = 60, height = 55, units = "mm", bg = "white"
)
ggsave(
  file.path("results/figures", "Figure1B.pdf"),
  plot = bplot, device = "pdf",
  width = 60, height = 55, units = "mm", bg = "white"
)
ggsave(
  file.path("results/figures", "Figure1C.pdf"),
  plot = cplot, device = "pdf",
  width = 60, height = 55, units = "mm", bg = "white"
)
