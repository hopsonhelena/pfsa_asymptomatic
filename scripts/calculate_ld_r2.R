# Calculate linkage disequilibrium from unmixed parasite genotypes.
calculate_ld_r2 <- function(data, group, pfsa1_column, pfsa3_column) {
  unmixed <- data[
    data$Infected == 1 &
      data[[pfsa1_column]] %in% c("0/0", "1/1") &
      data[[pfsa3_column]] %in% c("0/0", "1/1"),
    ,
    drop = FALSE
  ]

  pfsa1_alt <- as.integer(unmixed[[pfsa1_column]] == "1/1")
  pfsa3_alt <- as.integer(unmixed[[pfsa3_column]] == "1/1")
  pfsa1_alt_frequency <- mean(pfsa1_alt)
  pfsa3_alt_frequency <- mean(pfsa3_alt)
  joint_alt_frequency <- mean(pfsa1_alt == 1 & pfsa3_alt == 1)
  d <- joint_alt_frequency - pfsa1_alt_frequency * pfsa3_alt_frequency
  d_max <- min(
    pfsa1_alt_frequency * (1 - pfsa3_alt_frequency),
    pfsa3_alt_frequency * (1 - pfsa1_alt_frequency)
  )
  pearson_r <- cor(pfsa1_alt, pfsa3_alt, use = "pairwise.complete.obs")

  data.frame(
    Group = group,
    N = nrow(unmixed),
    Pfsa1_alt_frequency = pfsa1_alt_frequency,
    Pfsa3_alt_frequency = pfsa3_alt_frequency,
    Joint_alt_frequency = joint_alt_frequency,
    D = d,
    D_max = d_max,
    D_prime = d / d_max,
    r_squared = pearson_r^2
  )
}
