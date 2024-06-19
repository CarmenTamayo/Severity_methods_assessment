# Function to process simulation results
tidy_sim_results <- function(sim_data) {
  # Converting class to character so that it's possible to use ifelse
  sim_data$date_outcome <- as.character(sim_data$date_outcome)
  # Loop through each individual and update columns
  for (i in 1:nrow(sim_data)) {
    sim_data$date_death[i] <- ifelse(sim_data$outcome[i] == "died", sim_data$date_outcome[i], NA)
    sim_data$date_recovery[i] <- ifelse(sim_data$outcome[i] == "recovered", sim_data$date_outcome[i], NA)
  }
  # Convert date columns to Date class
  sim_data$date_death <- as.Date(sim_data$date_death, format = "%Y-%m-%d")
  sim_data$date_recovery <- as.Date(sim_data$date_recovery, format = "%Y-%m-%d")

  return(sim_data)
}

