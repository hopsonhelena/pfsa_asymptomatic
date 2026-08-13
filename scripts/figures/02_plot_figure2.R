# Plot main Figure 2 from its analysis tables.

library(ggplot2)
library(dplyr)
library(patchwork)

fig2ab <- read.delim("results/figures/data/Figure2AB_genotype_counts.tsv")
fig2cd <- read.delim("results/figures/data/Figure2CD_VAF.tsv")

theme_set(theme_minimal(base_size = 7))

counts <- fig2ab
vaf <- fig2cd
counts$HbS_genotype <- factor(counts$HbS_genotype, c("AA", "AS"))
counts$Parasite_genotype <- factor(
  counts$Parasite_genotype, c("0/0", "0/1", "1/1")
)

make_genotype_plot <- function(site) {
  plot_data <- counts[counts$Site == site, , drop = FALSE]
  ggplot(
    plot_data,
    aes(
      x = HbS_genotype, y = proportion, fill = Parasite_genotype
    )
  ) +
    geom_col(width = 0.6, color = "black", linewidth = 0.2) +
    scale_fill_manual(
      values = c(
        "0/0" = "#FEE0D2", "0/1" = "#FC9272", "1/1" = "#CB181D"
      ),
      labels = c(
        paste0(site, "-"), paste0(site, "+/-"), paste0(site, "+")
      ),
      name = NULL
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(y = "Proportion", x = NULL) +
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
      legend.title = element_text(size = 7)
    )
}

dplot <- make_genotype_plot("Pfsa1")
eplot <- make_genotype_plot("Pfsa3")

make_vaf_plot <- function(site, show_y = TRUE, show_legend = TRUE) {
  plot_data <- vaf[vaf$Site == site, , drop = FALSE]
  plot_data$HbS_genotype <- factor(
    plot_data$HbS_genotype, c("AA", "AS")
  )
  plot <- ggplot(
    plot_data,
    aes(x = VAF, color = HbS_genotype, fill = HbS_genotype)
  ) +
    geom_density(alpha = 0.5) +
    labs(
      x = paste0(site, " VAF"),
      color = NULL,
      fill = NULL,
      y = if (show_y) "Density" else NULL
    ) +
    theme_minimal(base_size = 7) +
    scale_color_manual(
      values = c("AA" = "deepskyblue4", "AS" = "lightgoldenrod2")
    ) +
    scale_fill_manual(
      values = c("AA" = "deepskyblue4", "AS" = "lightgoldenrod2")
    ) +
    theme(
      text = element_text(size = 7),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      legend.text = element_text(size = 7),
      legend.key.size = grid::unit(0.3, "cm")
    ) +
    coord_cartesian(ylim = c(0, 2.8))
  if (!show_legend) {
    plot <- plot + theme(legend.position = "none")
  }
  plot
}

fplot <- make_vaf_plot("Pfsa1", show_y = TRUE, show_legend = FALSE)
gplot <- make_vaf_plot("Pfsa3", show_y = FALSE, show_legend = TRUE)

top <- (dplot + eplot) + plot_layout(axis_titles = "collect")
bottom <- (fplot + gplot) +
  plot_layout(axis_titles = "collect", guides = "collect")
figure2 <- top / bottom + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 7))

dir.create("results/figures/plots", recursive = TRUE, showWarnings = FALSE)
ggsave(
  file.path("results/figures/plots", "Figure2.pdf"),
  plot = figure2, device = "pdf",
  width = 178, height = 127, units = "mm", bg = "white"
)
ggsave(
  file.path("results/figures/plots", "Figure2A.pdf"),
  plot = dplot, device = "pdf",
  width = 85, height = 65, units = "mm", bg = "white"
)
ggsave(
  file.path("results/figures/plots", "Figure2B.pdf"),
  plot = eplot, device = "pdf",
  width = 85, height = 65, units = "mm", bg = "white"
)
ggsave(
  file.path("results/figures/plots", "Figure2C.pdf"),
  plot = fplot, device = "pdf",
  width = 85, height = 65, units = "mm", bg = "white"
)
ggsave(
  file.path("results/figures/plots", "Figure2D.pdf"),
  plot = gplot, device = "pdf",
  width = 85, height = 65, units = "mm", bg = "white"
)
