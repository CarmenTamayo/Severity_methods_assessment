# Function to convert to incidence
convert_to_incidence <- function(sim) {
  incidence_data <- incidence2::incidence(sim,
                                          date_index = c("date_onset", "date_death", "date_recovery"),
                                          interval = 1L) %>%
    complete_dates() %>%
    pivot_wider(names_from = count_variable, values_from = count)
  #Vector with names for renaming
  rename_vec <- c("date_index" = "date", "date_onset" = "cases", "date_death" = "deaths" , "date_recovery" = "recovered")
  # Rename columns conditionally
  for (old_name in names(rename_vec)) {
    if (old_name %in% colnames(incidence_data)) {
      colnames(incidence_data)[colnames(incidence_data) == old_name] <- rename_vec[[old_name]]
    }
  }
  return(incidence_data)
}
