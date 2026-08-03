# Create and plot data for main Figure 4.

library(cowplot)
library(ggplot2)


# ---- Create table for Figure 4 ----

metadata <- read.table("inputs/metadata.tsv", check.names = FALSE, header = TRUE, sep = "\t")
pfsa1_column <- "Pf3D7_02_v3.631190_gt"
pfsa3_column <- "Pf3D7_11_v3.1058035_gt"

unmixed <- metadata[
  metadata$Infected == 1 &
    metadata[[pfsa1_column]] %in% c("0/0", "1/1") &
    metadata[[pfsa3_column]] %in% c("0/0", "1/1"),
  ,
  drop = FALSE
]
pfsa1_alt_frequency <- mean(unmixed[[pfsa1_column]] == "1/1")
pfsa3_alt_frequency <- mean(unmixed[[pfsa3_column]] == "1/1")
expected_probability <- c(
  (1 - pfsa1_alt_frequency) * (1 - pfsa3_alt_frequency),
  (1 - pfsa1_alt_frequency) * pfsa3_alt_frequency,
  pfsa1_alt_frequency * (1 - pfsa3_alt_frequency),
  pfsa1_alt_frequency * pfsa3_alt_frequency
)
combinations <- data.frame(
  Pfsa1_genotype = c("0/0", "0/0", "1/1", "1/1"),
  Pfsa3_genotype = c("0/0", "1/1", "0/0", "1/1")
)
get_joint_counts <- function(data, group) {
  observed <- vapply(seq_len(nrow(combinations)), function(i) {
    sum(
      data[[pfsa1_column]] == combinations$Pfsa1_genotype[i] &
        data[[pfsa3_column]] == combinations$Pfsa3_genotype[i]
    )
  }, numeric(1))
  data.frame(
    Group = group,
    N = nrow(data),
    combinations,
    Observed_count = observed,
    Expected_count = nrow(data) * expected_probability,
    Observed_proportion = observed / nrow(data),
    Expected_proportion = expected_probability,
    Pfsa1_alt_frequency_for_null = pfsa1_alt_frequency,
    Pfsa3_alt_frequency_for_null = pfsa3_alt_frequency
  )
}
fig4a <- rbind(
  get_joint_counts(unmixed, "All"),
  get_joint_counts(unmixed[unmixed$HbS_gt == "AA", ], "AA"),
  get_joint_counts(unmixed[unmixed$HbS_gt == "AS", ], "AS")
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
fig4b <- mixed_either[, c(
  "Seq.Sample.ID", "Sample.ID", "HbS_gt", "pfsa1_vaf", "pfsa3_vaf"
)]
names(fig4b) <- c(
  "Sequencing_sample_ID", "Study_sample_ID", "HbS_genotype",
  "Pfsa1_VAF", "Pfsa3_VAF"
)

write.table(
  fig4a, "figure_data/Figure4A_observed_expected.tsv", sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  fig4b, "figure_data/Figure4B_paired_VAF.tsv", sep = "\t", row.names = FALSE, quote = FALSE
)

# ---- Plot figure ----

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
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

ggsave(
  file.path("figures", "Figure4A.pdf"),
  plot = ld_row, device = "pdf",
  width = 190, height = 105, units = "mm", bg = "white"
)
ggsave(
  file.path("figures", "Figure4B.pdf"),
  plot = vaf_row, device = "pdf",
  width = 150, height = 75, units = "mm", bg = "white"
)
