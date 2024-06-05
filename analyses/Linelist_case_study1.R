
# Script for individual-level data case study 1: Linelist with known outcomes,
# with a similar delay between onset to death and onset to recovery (this is
# relevant because the onset-death delay will be used to estimate the % of known
# outcomes). 3 different approaches for estimating CFR in real-time during an
# outbreak will be used and compared to the "true" CFR estimate at the end of
# the outbreak (i.e., when all cases have had a registered outcome):
# 1) Naive estimate of deaths/cases
# 2) Retrospective cohort by omitting latest cases for which we don't expect to have an outcome (i.e., within 1 onset-death before real-time)
# 3) Prospective cohort by taking into account deaths corresponding to cases in rt (i.e., within 1 onset-death after real-time)
# 4) Using the {cfr} function `known_outcomes()` to estimate the no. of cases
# with a known outcome and use this as the denominator.

# Using simulist to simulate outbreak data
# Columns that we need: onset date, column to indicate outcome, outcome date

library(simulist)
library(epiparameter)
library(tidyverse)
library(incidence2)
library(cfr)
library(data.table)

contact_distribution <- epidist(
  disease = "MERS",
  epi_dist = "contact distribution",
  prob_distribution = "pois",
  prob_distribution_params = c(mean = 3)
)

ip_mers <- epidist(
  disease = "MERS",
  epi_dist = "infectious_period",
  pathogen = "MERS-CoV",
  prob_distribution = "gamma",
  prob_distribution_params = c(shape = 3, scale = 3)
)

o_d_mers <- epidist_db(
  disease = "MERS",
  epi_dist = "onset to death",
  single_epidist = TRUE
)

o_r_mers <- onset_to_death <- epidist_db(
  disease = "MERS",
  epi_dist = "onset to death",
  single_epidist = TRUE
)

#Note: does it make sense to use a distribution from epiparameter to simulate the linelist
#and then use this same delay in known_outcomes? that's not very realistic

#Can we choose the length of the outbreak as well as the size?

set.seed(12345678)
linelist <- sim_linelist(
  contact_distribution = contact_distribution,
  infect_period = ip_mers,
  onset_to_recovery = o_r_mers,
  prob_infect = 0.7,
  onset_to_hosp = NA,
  hosp_risk = NA,
  non_hosp_death_risk = 0.3,
  onset_to_death = o_d_mers,
  outbreak_size = c(500, 1000)
)


#Creating columns with dates of death and recovery
linelist$date_death <- fifelse(linelist$outcome == "died", linelist$date_outcome, NA)
linelist$date_recovery <- fifelse(linelist$outcome == "recovered", linelist$date_outcome, NA)

linelist_mers <- linelist[,c(1,6,13,14)]

#Visualising to choose a cut-off point for the outbreak
incidence_mers <- incidence2::incidence(linelist_mers, date_index = c("date_onset", "date_death", "date_recovery"), interval = 1L)
incidence_recovered <- incidence_mers %>% filter(count_variable == "date_recovery") %>% select(date = date_index, recovered = count)
incidence_mers_long <- cfr::prepare_data(incidence_mers, cases_variable = "date_onset", deaths_variable = "date_death")

incidence_mers_long <- incidence_mers_long %>%
  left_join(incidence_recovered, by= "date") %>%
  mutate(recovered = replace_na(recovered, 0))


incidence_mers_long$date <- as.Date(incidence_mers_long$date)
ggplot(incidence_mers_long, aes(x = date)) + geom_point(aes(y = cases), colour = "blue") +
  geom_point(aes(y = deaths), colour = "red") + geom_point(aes(y = recovered), colour = "green") + theme_bw() + scale_x_date(date_breaks = "7 days")

#Cut-off point for real-time outbreak #1- During outbreak peak

real_time <- "2023-02-09"
outbreak_mers_rt <- incidence_mers_long[incidence_mers_long$date <= real_time,]

#"True" CFR, estimated at the end of the outbreak as total deaths/total cases

true_CFR <- sum(incidence_mers_long$deaths)/sum(incidence_mers_long$cases)

#Option 1- Naive CFR in real-time

naive_CFR <- sum(outbreak_mers_rt$deaths)/sum(outbreak_mers_rt$cases)

#Option 2- Retrospective cohort where onsets closest to current time are removed

onset_to_death$summary_stat$mean
#Taking into account that the mean onset-death is around 2 weeks, we only keep those
#cases whose onset is up to 2 weeks before

cutoff_1 <- max(outbreak_mers_rt$date) - 14
outbreak_mers_cohort1 <- outbreak_mers_rt

outbreak_mers_cohort1$cases <- ifelse(outbreak_mers_cohort1$date > cutoff_1, 0, outbreak_mers_cohort1$cases)
retrospective_cfr <- sum(outbreak_mers_cohort1$deaths)/sum(outbreak_mers_cohort1$cases)

#Option 3- Forward-looking cohort that includes deaths that occur 1 onset-death after current time onsets

cutoff_deaths1 <- max(outbreak_mers_rt$date) + 14
outbreak_mers_cohort_deaths1 <- incidence_mers_long[incidence_mers_long$date <= cutoff_deaths1,]

forward_cfr <- sum(outbreak_mers_cohort_deaths1$deaths)/sum(outbreak_mers_rt$cases)

#Option 4- Using onset-death delay to estimate known outcomes

od_parameters <- get_parameters(o_d_mers)
known_outcomes <- cfr::estimate_outcomes(data = outbreak_mers_rt, delay_density = function(x) dgamma(x, shape = od_parameters[["shape"]], scale = od_parameters[["scale"]]))
cfr_known_outcomes <- sum(outbreak_mers_rt$deaths)/sum(known_outcomes$estimated_outcomes)


# On the growing phase of the epidemic, option 1 severely underestimates CFR, option 2 overestimates CFR,
# option 3 provides an overestimated CFR, and option 4 provides an almost exact CFR estimate













