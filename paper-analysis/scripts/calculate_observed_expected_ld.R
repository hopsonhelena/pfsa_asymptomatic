# Calculate observed and expected counts and proportions, as well as binomial p-values for each Pfsa1 and Pfsa3 genotype combination in unmixed infections.
calculate_observed_expected_ld <- function(
  data, all_sample_data, group, pfsa1_column, pfsa3_column
) {
  unmixed <- data[
    data$Infected == 1 &
      data[[pfsa1_column]] %in% c("0/0", "1/1") &
      data[[pfsa3_column]] %in% c("0/0", "1/1"),
    ,
    drop = FALSE
  ]

  all_samples_unmixed <- all_sample_data[
    all_sample_data$Infected == 1 &
      all_sample_data[[pfsa1_column]] %in% c("0/0", "1/1") &
      all_sample_data[[pfsa3_column]] %in% c("0/0", "1/1"),
    ,
    drop = FALSE
  ]
  pfsa1_alt_frequency <- mean(
    all_samples_unmixed[[pfsa1_column]] == "1/1"
  )
  pfsa3_alt_frequency <- mean(
    all_samples_unmixed[[pfsa3_column]] == "1/1"
  )
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

  observed <- vapply(seq_len(nrow(combinations)), function(i) {
    sum(
      unmixed[[pfsa1_column]] == combinations$Pfsa1_genotype[i] &
        unmixed[[pfsa3_column]] == combinations$Pfsa3_genotype[i]
    )
  }, numeric(1))
  binomial_p_value <- vapply(seq_len(nrow(combinations)), function(i) {
    binom.test(
      observed[i], nrow(unmixed), p = expected_probability[i]
    )$p.value
  }, numeric(1))
  data.frame(
    Group = group,
    N = nrow(unmixed),
    combinations,
    Observed_count = observed,
    Expected_count = nrow(unmixed) * expected_probability,
    Binomial_p_value = binomial_p_value,
    Observed_proportion = observed / nrow(unmixed),
    Expected_proportion = expected_probability,
    Pfsa1_alt_frequency_for_null = pfsa1_alt_frequency,
    Pfsa3_alt_frequency_for_null = pfsa3_alt_frequency
  )
}
