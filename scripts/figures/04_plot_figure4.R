# Plot main Figure 4 from its analysis tables.

library(cowplot)
library(ggplot2)

fig4a <- read.delim("results/figures/data/Figure4A_observed_expected.tsv")
fig4b <- read.delim("results/figures/data/Figure4B_paired_VAF.tsv")

theme_set(theme_minimal(base_size = 7))

make_ld_plot <- function(results, title) {
  x_sign <- c(-1, -1, 1, 1)
  y_sign <- c(-1, 1, -1, 1)
  square_data <- data.frame(
    xmin = 0,
    xmax = rep(x_sign, 2) *
      sqrt(c(results$Observed_count, results$Expected_count)),
    ymin = 0,
    ymax = rep(y_sign, 2) *
      sqrt(c(results$Observed_count, results$Expected_count)),
    fill = c(rep("grey75", 4), rep(NA_character_, 4)),
    line_color = c(rep("black", 4), rep("blue", 4)),
    line_type = c(rep("solid", 4), rep("dashed", 4)),
    stringsAsFactors = FALSE
  )

  labels <- paste0(
    "Pfsa1", c("-", "-", "+", "+"),
    "Pfsa3", c("-", "+", "-", "+"),
    "\nn=", results$Observed_count
  )

  ggplot(square_data) +
    geom_rect(
      aes(
        xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = fill, color = line_color, linetype = line_type
      )
    ) +
    scale_fill_identity(na.value = "transparent") +
    scale_color_identity() +
    scale_linetype_identity() +
    geom_vline(xintercept = 0, color = "black") +
    geom_hline(yintercept = 0, color = "black") +
    annotate(
      "text",
      x = c(-14, -14, 14, 14),
      y = c(-27, 27, -27, 27),
      label = labels,
      size = 7 / .pt
      ) +
    coord_fixed(xlim = c(-28, 28), ylim = c(-29, 28)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 7) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 7)
    )
}

all_plot <- make_ld_plot(fig4a[fig4a$Group == "All", ], "All")
aa_plot <- make_ld_plot(fig4a[fig4a$Group == "AA", ], "AA")
as_plot <- make_ld_plot(fig4a[fig4a$Group == "AS", ], "AS")

fig4b$HbS_genotype <- factor(fig4b$HbS_genotype, c("AA", "AS"))
make_vaf_plot <- function(data, genotype, show_y_title = TRUE) {
  plot_data <- data[data$HbS_genotype == genotype, , drop = FALSE]
  ggplot(
    plot_data,
    aes(x = Pfsa1_VAF, y = Pfsa3_VAF)
  ) +
    geom_point(size = 1, color = "red", shape = 15) +
    labs(
      title = genotype,
      x = "Pfsa1 VAF",
      y = if (show_y_title) "Pfsa3 VAF" else NULL
    ) +
    theme_minimal(base_size = 7) +
    theme(legend.position = "none")
}
vaf_plot_aa <- make_vaf_plot(fig4b, "AA")
vaf_plot_as <- make_vaf_plot(fig4b, "AS", show_y_title = FALSE)

ld_row <- plot_grid(
  all_plot, aa_plot, as_plot,
  nrow = 1,
  align = "h"
)
vaf_row <- plot_grid(
  vaf_plot_aa, vaf_plot_as,
  nrow = 1,
  align = "h",
  axis = "tb"
)

figure4_body <- plot_grid(
  ld_row, vaf_row,
  ncol = 1,
  rel_heights = c(0.8, 0.5)
)
figure4 <- ggdraw(figure4_body) +
  draw_plot_label(
    label = c("a", "b"),
    x = c(0.01, 0.01),
    y = c(0.99, 0.45),
    fontface = "bold",
    size = 7
  )

dir.create("results/figures/plots", recursive = TRUE, showWarnings = FALSE)
ggsave(
  file.path("results/figures/plots", "Figure4.pdf"),
  plot = figure4, device = "pdf",
  width = 180, height = 168.6, units = "mm", bg = "white"
)
ggsave(
  file.path("results/figures/plots", "Figure4A.pdf"),
  plot = ld_row, device = "pdf",
  width = 180, height = 99.5, units = "mm", bg = "white"
)
ggsave(
  file.path("results/figures/plots", "Figure4B.pdf"),
  plot = vaf_row, device = "pdf",
  width = 150, height = 75, units = "mm", bg = "white"
)
