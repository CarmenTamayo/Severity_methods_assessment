# Function to select a random date within an outbreak phase and create a subset
# of data between relevant dates of the outbreak
subset_outbreak_phase <- function(data, from, to) {
  # Sample a random date between "from" and "to" statuses
  sampled_date <- data %>%
    slice(sample(seq(which(outbreak_dates == from), which(outbreak_dates == to)), 1)) %>%
    pull(date) %>%
    as.character()

  # Subset data based on the sampled date
  subset_data <- data[data$date <= sampled_date, ]

  return(subset_data)
}
