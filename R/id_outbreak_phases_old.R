#Function to smooth outbreak trend and create a column that indicates the outbreak's
# peak (mean of window with highest no. of cases), start, and end
id_outbreak_phases_old <- function(inc_data) {
  # setting window size
  window_size <- 7
  trend <- inc_data %>%
    mutate(cases_ma = rollmean(cases, k = window_size, fill = NA, align = "center")) %>%
    mutate(date = as.Date(date)) %>%
    mutate(cases_ma = replace_na(cases_ma, 0)) %>%
    mutate(outbreak_dates = case_when(
      row_number() == which.max(cases_ma) ~ "peak",
      row_number() == which.min(if_else(cases_ma > 1, row_number(), Inf)) ~ "start",
      date == max(date) ~ "end",
      TRUE ~ NA_character_
    ))
}


