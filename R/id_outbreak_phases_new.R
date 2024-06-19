# Function that takes the outbreak's peak (identified previously by id_outbreak_peak),
# pastes this on a new column

id_outbreak_phases_new <- function(inc_data, outbreak_peaks){
  inc_data$date <- as.Date(inc_data$date)
  # Id start date
  start_date <- min(inc_data$date)
  # Select a possible date for the growing phase
  possible_growing_dates <- inc_data %>%
    filter(date > start_date + 10 & date < outbreak_peaks) %>%
    pull(date)
  # Select one of these dates at random
  growing_date <- sample(possible_growing_dates, 1)
  # Add new column with outbreak dates
  inc_data <- inc_data %>%
    mutate(outbreak_dates = case_when(
      date == outbreak_peaks & !duplicated(date) ~ "peak",
      cases >= 5 & row_number() == min(which(cases >= 5))  ~ "start",
      date == max(date) & !duplicated(date) ~ "end",
      date == growing_date ~ "growing",
      TRUE ~ NA_character_
    )
    )
  return(inc_data)
}
