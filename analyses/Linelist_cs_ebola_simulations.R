
# Script for individual-level data case study 4: disease with a delay from onset
# to death shorter than onset to recovery (Ebola-like outbreak).
# Methods to compare:
# 1) Naive estimate of deaths/cases
# 2) Retrospective cohort by omitting latest cases for which we don't expect to have an outcome (i.e., within 1 onset-death before real-time)
# 3) Using {cfr} function estimate_outcomes to recalculate denominator
# 4) Using {cfr} function cfr_static


# Using simulist to simulate outbreak data
# Columns that we need: onset date, column to indicate outcome, outcome date (recovery/death)

library(data.table)
library(zoo)
library(tidyverse)
library(simulist)
library(epiparameter)
library(incidence2)
library(kableExtra)

contact_distribution <- epidist(
  disease = "ebola",
  epi_dist = "contact distribution",
  prob_distribution = "pois",
  prob_distribution_params = c(mean = 3)
)

ebola_ip_meansd <- c(mean = 9.4, sd = 5.5)
# From https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4560911/

ebola_ip_params <- convert_summary_stats_to_params("lnorm", mean = ebola_ip_meansd[["mean"]], sd = ebola_ip_meansd[["sd"]])

ip_ebola_epidist <- epidist(
  disease = "ebola",
  epi_dist = "infectious_period",
  prob_distribution = "lnorm",
  prob_distribution_params = c(meanlog = ebola_ip_params$meanlog, sdlog = ebola_ip_params$sdlog),
  auto_calc_params = TRUE
)

# No distributions in epireview
# summary statistics: mean = 18 and sd = 4

onset_recovery_ebola <- convert_summary_stats_to_params("gamma", mean = 18, sd = 4)

o_r_ebola_epidist <- epidist(
  disease = "ebola",
  epi_dist = "onset to recovery",
  prob_distribution = "gamma",
  prob_distribution_params = c(shape = onset_recovery_ebola[["shape"]], scale = onset_recovery_ebola[["scale"]])
)

o_d_ebola_epidist <- epidist_db(
  disease = "ebola",
  epi_dist = "onset to death",
  single_epidist = TRUE
)

# Function rum_ebola_simulations uses simulist to create a list with n=100 (for now) linelists
set.seed(487593478)
simulations_ebola <- run_ebola_simulations(100)

# Function tidy_sim_results creates 2 new columns (date_onset and date_recovery) and makes sure they're on a date format
linelist_with_added_columns <- lapply(simulations_ebola, tidy_sim_results)


# Function convert_to_incidence uses incidence2 to obtain daily incidence data and modify format so it matches {cfr} requirements
incidence_ebola <- lapply(linelist_with_added_columns, convert_to_incidence)


# Identifying outbreak peak
outbreak_peaks <- sapply(linelist_with_added_columns, id_outbreak_peak)

# Function id_outbreak_phases to establish the trend of the outbreak and identify the start, peak, and end-points

incidence_ebola_trend <- Map(id_outbreak_phases_new, incidence_ebola, outbreak_peaks)

# Function subset_outbreak_phase to create 2 lists of dfs, one with a random cut-off point in the growing phase,
# and another one with a random cut-off point in the declining phase of the outbreak
rt_ebola_growing <- lapply(incidence_ebola_trend, subset_outbreak_phase, from = "start", to = "peak")
rt_ebola_decline <- lapply(incidence_ebola_trend, subset_outbreak_phase, from = "peak", to = "end")

# 3rd list of dfs, with the cut-off point at the peak of the outbreak
# Ideally this would be done with subset_outbreak_phase as well
rt_ebola_peak <- incidence_ebola_trend
for (i in seq_along(rt_ebola_peak)) {
  peak_index <- which(rt_ebola_peak[[i]]$outbreak_dates == "peak")
  rt_ebola_peak[[i]] <- rt_ebola_peak[[i]][1:peak_index, ]
}

# True CFR at the end of the outbreak:

true_cfrs <- sapply(incidence_ebola, estimate_cfr, numerator = "deaths", denominator = "cases")
true_cfr_ci <- cfr_ci(true_cfrs)

### Naive CFRs
# Growing phase
naive_cfr_growing <- sapply(rt_ebola_growing, estimate_cfr, numerator = "deaths", denominator = "cases")
naive_cfr_growing_ci <- cfr_ci(naive_cfr_growing)

# Peak
naive_cfr_peak <- sapply(rt_ebola_peak, estimate_cfr, numerator = "deaths", denominator = "cases")
naive_cfr_peak_ci <- cfr_ci(naive_cfr_peak)

# Decline
naive_cfr_decline <- sapply(rt_ebola_decline, estimate_cfr, numerator = "deaths", denominator = "cases")
naive_cfr_decline_ci <- cfr_ci(naive_cfr_decline)

### Retrospective cohort
# Substitute the column "cases" by 0 if the column "date" is within 8 days of the last/maximum date

# Growing phase
# Here taking the growing phase cutoff as a random date between the row that says
# "growing", which is at least 1 o-d away from the start of the outbreak, to
# "peak"- "start" cannot be used since there needs to be enough deaths/cases
# to create a retrospective cohort

rt_ebola_growing_cohort <- lapply(incidence_ebola_trend, subset_outbreak_phase, from = "growing", to = "peak")

for (i in seq_along(rt_ebola_growing_cohort)) {
  df <- rt_ebola_growing_cohort[[i]]
  df$cases <- ifelse(max(df$date)-df$date <= 8, 0, df$cases)
  rt_ebola_growing_cohort[[i]] <- df
}

rc_growing_cfr <- sapply(rt_ebola_growing_cohort, estimate_cfr, numerator = "deaths", denominator = "cases")
rc_growing_cfr_ci <- cfr_ci(rc_growing_cfr)


# Peak

rt_ebola_peak_cohort <- rt_ebola_peak

for (i in seq_along(rt_ebola_peak_cohort)) {
  df <- rt_ebola_peak_cohort[[i]]
  df$cases <- ifelse(max(df$date)-df$date <= 8, 0, df$cases)
  rt_ebola_peak_cohort[[i]] <- df
}

rc_peak_cfr <- sapply(rt_ebola_peak_cohort, estimate_cfr, numerator = "deaths", denominator = "cases")
rc_peak_ci <- cfr_ci(rc_peak_cfr)


# Decline

rt_ebola_decline_cohort <- rt_ebola_decline

for (i in seq_along(rt_ebola_decline_cohort)) {
  df <- rt_ebola_decline_cohort[[i]]
  df$cases <- ifelse(max(df$date)-df$date <= 8, 0, df$cases)
  rt_ebola_decline_cohort[[i]] <- df
}

rc_decline_cfr <- sapply(rt_ebola_decline_cohort, estimate_cfr, numerator = "deaths", denominator = "cases")
rc_decline_ci <- cfr_ci(rc_decline_cfr)

## Delay-adjusted ratio
# Getting onset-death delay from {epiparameter}

o_d_ebola_epidist_extracted <- epidist_db(
  disease = "ebola",
  epi_dist = "onset to death",
  author = "The Ebola Outbreak Epidemiology Team"
)
o_d_ebola_params <- get_parameters(o_d_ebola_epidist_extracted)

# Growing phase- here we can take any date from the beginning of the outbreak to the peak
# because the delay-adjusted method can be used with any no. deaths, even when there
# haven't been any deaths yet

known_outcomes_growing <- lapply(rt_ebola_growing, cfr::estimate_outcomes, function(x) dgamma(x, shape = o_d_ebola_params[["shape"]], scale = o_d_ebola_params[["scale"]]))
cfr_delay_growing <- sapply(known_outcomes_growing, estimate_cfr, numerator = "deaths", denominator = "estimated_outcomes")
cfr_delay_growing_ci <- cfr_ci(cfr_delay_growing)

# Peak

known_outcomes_peak <- lapply(rt_ebola_peak, cfr::estimate_outcomes, function(x) dgamma(x, shape = o_d_ebola_params[["shape"]], scale = o_d_ebola_params[["scale"]]))
cfr_delay_peak <- sapply(known_outcomes_peak, estimate_cfr, numerator = "deaths", denominator = "estimated_outcomes")
cfr_delay_peak_ci <- cfr_ci(cfr_delay_peak)

# Decline

known_outcomes_decline <- lapply(rt_ebola_decline, cfr::estimate_outcomes, function(x) dgamma(x, shape = o_d_ebola_params[["shape"]], scale = o_d_ebola_params[["scale"]]))
cfr_delay_decline <- sapply(known_outcomes_decline, estimate_cfr, numerator = "deaths", denominator = "estimated_outcomes")
cfr_delay_decline_ci <- cfr_ci(cfr_delay_decline)


## Using {cfr_static}

# Growing

cfr_cfr_growing <- unlist(sapply(rt_ebola_growing, function(x) {
  result <- cfr::cfr_static(x, function(x) dgamma(x, shape = o_d_ebola_params[["shape"]], scale = o_d_ebola_params[["scale"]]))
  result[, 1] * 100  # Extracting only the first column that contains the mean cfr
}))

cfr_cfr_growing_ci <- cfr_ci(cfr_cfr_growing)

# Peak

cfr_cfr_peak <- unlist(sapply(rt_ebola_peak, function(x) {
  result <- cfr::cfr_static(x, function(x) dgamma(x, shape = o_d_ebola_params[["shape"]], scale = o_d_ebola_params[["scale"]]))
  result[, 1] * 100  # Extracting only the first column that contains the mean cfr
}))

cfr_cfr_peak_ci <- cfr_ci(cfr_cfr_peak)

# Decline

cfr_cfr_decline <- unlist(sapply(rt_ebola_decline, function(x) {
  result <- cfr::cfr_static(x, function(x) dgamma(x, shape = o_d_ebola_params[["shape"]], scale = o_d_ebola_params[["scale"]]))
  result[, 1] * 100  # Extracting only the first column that contains the mean cfr
}))

cfr_cfr_decline_ci <- cfr_ci(cfr_cfr_decline)



# Making table

table_cs_4 <- data.frame(
  Method = c("Naive estimate", "Retrospective cohort", "Delay-adjusted ratio", "Using {cfr}", "True CFR"),
  Growing = c(naive_cfr_growing_ci$mean, rc_growing_cfr_ci$mean, cfr_delay_growing_ci$mean, cfr_cfr_growing_ci$mean, true_cfr_ci$mean),
  `95% CI` = c(
    paste(naive_cfr_growing_ci$lower_ci, naive_cfr_growing_ci$upper_ci, sep = ", "),
    paste(rc_growing_cfr_ci$lower_ci, rc_growing_cfr_ci$upper_ci, sep = ", "),
    paste(cfr_delay_growing_ci$lower_ci, cfr_delay_growing_ci$upper_ci, sep = ", "),
    paste(cfr_cfr_growing_ci$lower_ci, cfr_cfr_growing_ci$upper_ci, sep = ", "),
    paste(true_cfr_ci$lower_ci, true_cfr_ci$upper_ci, sep = ", ")),
  Peak = c(naive_cfr_peak_ci$mean, rc_peak_ci$mean, cfr_delay_peak_ci$mean, cfr_cfr_peak_ci$mean, true_cfr_ci$mean),
  `95% CI (Peak)` = c(
    paste(naive_cfr_peak_ci$lower_ci, naive_cfr_peak_ci$upper_ci, sep = ", "),
    paste(rc_peak_ci$lower_ci, rc_peak_ci$upper_ci, sep = ", "),
    paste(cfr_delay_peak_ci$lower_ci, cfr_delay_peak_ci$upper_ci, sep = ", "),
    paste(cfr_cfr_peak_ci$lower_ci, cfr_cfr_peak_ci$upper_ci, sep = ", "),
    paste(true_cfr_ci$lower_ci, true_cfr_ci$upper_ci, sep = ", ")
  ),
  Declining = c(naive_cfr_decline_ci$mean, rc_decline_ci$mean, cfr_delay_decline_ci$mean, cfr_cfr_decline_ci$mean, true_cfr_ci$mean),
  `95% CI (Declining)` = c(
    paste(naive_cfr_decline_ci$lower_ci, naive_cfr_decline_ci$upper_ci, sep = ", "),
    paste(rc_decline_ci$lower_ci, rc_decline_ci$upper_ci, sep = ", "),
    paste(cfr_delay_decline_ci$lower_ci, cfr_delay_decline_ci$upper_ci, sep = ", "),
    paste(cfr_cfr_decline_ci$lower_ci, cfr_cfr_decline_ci$upper_ci, sep = ", "),
    paste(true_cfr_ci$lower_ci, true_cfr_ci$upper_ci, sep = ", ")
  )
  )


color_func <- function(value, threshold) {
  if (abs(value - threshold) <= 0.001) {
    "white"
  } else if (abs(value - threshold) <= 2) {
    "lightgreen"
  } else if (value < threshold) {
    "lightblue"
  } else if (value > threshold) {
    "lightcoral"
  }
}

cols_to_color <- c("Growing", "Peak", "Declining")

cfr_table_colours <- table_cs_4
for (col in cols_to_color) {
  cfr_table_colours[[col]] <- sapply(cfr_table_colours[[col]], function(x) {
    paste0("<span style='background-color:", color_func(x, true_cfr_ci$mean), "; display: block;'>", x, "</span>")
  })
}

# Print the data frame with coloured cells
kable(cfr_table_colours, format = "html", escape = FALSE, col.names = c("Method", "Growing phase CFR mean", "95% CI", "Peak CFR mean", "95% CI", "Decline phase CFR mean", "95% CI")) %>%
  kable_styling(full_width = FALSE) %>%
  row_spec(5, bold = TRUE)




