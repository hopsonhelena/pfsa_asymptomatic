# Create and plot data for main Figure 2.

library(ggplot2)
library(dplyr)


# ---- Create table for Figure 2 ----

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

fig2ab <- rbind(
  get_pfsa_counts(metadata, "Pfsa1", "Pf3D7_02_v3.631190_gt"),
  get_pfsa_counts(metadata, "Pfsa3", "Pf3D7_11_v3.1058035_gt")
)

get_mixed_vaf <- function(
  site, genotype_column, vaf_column
) {
  keep <- metadata$HbS_gt %in% c("AA", "AS") &
    metadata[[genotype_column]] == "0/1" &
    !is.na(metadata[[vaf_column]])
  output <- metadata[
    keep,
    c("Seq.Sample.ID", "HbS_gt", vaf_column),
    drop = FALSE
  ]
  names(output) <- c(
    "Sequencing_sample_ID", "HbS_genotype", "VAF"
  )
  output$Site <- site
  output[, c(
    "Site", "Sequencing_sample_ID", "HbS_genotype", "VAF"
  )]
}

fig2cd <- rbind(
  get_mixed_vaf(
    "Pfsa1", "Pf3D7_02_v3.631190_gt", "pfsa1_vaf"
  ),
  get_mixed_vaf(
    "Pfsa3", "Pf3D7_11_v3.1058035_gt", "pfsa3_vaf"
  )
)

write.table(
  fig2ab, "figure_data/Figure2AB_genotype_counts.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  fig2cd, "figure_data/Figure2CD_VAF.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

# ---- Plot figure ----

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
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
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
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

ggsave(
  file.path("figures", "Figure2A.pdf"),
  plot = dplot, device = "pdf",
  width = 85, height = 65, units = "mm", bg = "white"
)
ggsave(
  file.path("figures", "Figure2B.pdf"),
  plot = eplot, device = "pdf",
  width = 85, height = 65, units = "mm", bg = "white"
)
ggsave(
  file.path("figures", "Figure2C.pdf"),
  plot = fplot, device = "pdf",
  width = 85, height = 65, units = "mm", bg = "white"
)
ggsave(
  file.path("figures", "Figure2D.pdf"),
  plot = gplot, device = "pdf",
  width = 85, height = 65, units = "mm", bg = "white"
)
