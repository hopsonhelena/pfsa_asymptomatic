# Plot main Figure 3 from its analysis tables.

library(ggplot2)
library(patchwork)
library(dplyr)

group_levels <- c(
  "Gambia (severe)", "Kenya (severe)", "Ghana (mild)",
  "Cameroon (asymptomatic)"
)
variant_levels <- c("Pfsa1", "Pfsa3")

fig3a_output <- read.delim("results/tables/Figure3A_odds_ratios.tsv")
fig3b <- read.delim("results/tables/Figure3B_genotype_counts.tsv")

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

figure3_body <- top_row / bottom_row
figure3_body <- figure3_body &
  theme(plot.margin = margin(1, 0, 0, -20))
figure3 <- cowplot::ggdraw(figure3_body) +
  cowplot::draw_plot_label(
    label = c("a", "b"),
    x = c(0.05, 0.05),
    y = c(1.01, 0.5),
    fontface = "bold",
    size = 7
  )

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave(
  file.path("results/figures", "Figure3.pdf"),
  plot = figure3, device = "pdf",
  width = 178, height = 102, units = "mm", bg = "white"
)
ggsave(
  file.path("results/figures", "Figure3A.pdf"),
  plot = top_row, device = "pdf",
  width = 178, height = 50, units = "mm", bg = "white"
)
ggsave(
  file.path("results/figures", "Figure3B.pdf"),
  plot = bottom_row, device = "pdf",
  width = 178, height = 55, units = "mm", bg = "white"
)
