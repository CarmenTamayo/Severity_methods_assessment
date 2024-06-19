# Function to run multiple instances of sim_linelist
run_ebola_simulations <- function(sim_count) {
  # Initialize empty list to store simulation results
  results <- vector("list", length = sim_count)

  # Loop through simulations
  for (i in 1:sim_count) {
    # Set seed for reproducibility
    set.seed(i)  # Use loop counter for unique seeds

    # Run the simulation code and store the result
    results[[i]] <- sim_linelist(
      contact_distribution = contact_distribution,
      infectious_period = ip_ebola_epidist,
      onset_to_hosp = NA,
      hosp_risk = NA,
      hosp_death_risk = NA,
      onset_to_death = o_d_ebola_epidist,
      onset_to_recovery = o_r_ebola_epidist,
      prob_infect = 0.7,
      outbreak_size = c(500, 2000),
      non_hosp_death_risk = 0.3
    )
  }
  # Return the list of simulation results
  return(results)
}
