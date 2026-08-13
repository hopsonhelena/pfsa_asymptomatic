# Summarize a binary outcome within groups using a Wald 95% confidence interval.
count_rate_rows <- function(data, group_column, outcome_column) {
  groups <- split(data, data[[group_column]], drop = TRUE)
  output <- lapply(names(groups), function(group) {
    values <- groups[[group]][[outcome_column]]
    values <- values[!is.na(values)]
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
