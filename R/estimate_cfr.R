# Function to calculate cfr by specifying columns from a dataframe with
# deaths and cases
estimate_cfr <- function(df, numerator, denominator) {
  total_deaths <- sum(df[[numerator]])
  total_cases <- sum(df[[denominator]])
  cfr <- total_deaths/total_cases*100
  return(cfr)
}
