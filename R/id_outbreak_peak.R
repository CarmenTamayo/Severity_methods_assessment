# Function to estimate the peak of the outbreak using incidence2
# Incidence object from linelist has 3 count_variables, but we want the column
# date_onset only to be considered to estimate the peak
id_outbreak_peak <- function(linelist) {
  inc_data <- linelist %>%
    incidence2::incidence(date_index = "date_onset", interval = 1L) %>%
    complete_dates()
  peak <- inc_data %>%
    estimate_peak(first_only = TRUE) %>%
    pull(observed_peak) %>%
    as.character()
  return(peak)

}
