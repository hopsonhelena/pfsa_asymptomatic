# Count parasite genotypes and calculate proportions within HbS groups.
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
