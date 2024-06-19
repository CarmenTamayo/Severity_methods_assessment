# Function to estimate mean CFR with 95% CI from a numeric vector
cfr_ci <- function(x) {
  mean <- mean(x)
  sd <- sd(x)
  sem <- sd/ sqrt(length(x))
  margin_of_error <- sem * 1.96
  lower <- round(mean - margin_of_error, 2)
  upper <- round(mean + margin_of_error, 2)
  result <- list(
    mean = round(mean, 2),
    lower_ci = lower,
    upper_ci = upper
  )
  return(result)
}
