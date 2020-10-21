#------------------------------------------------------
# Calculate performance measures
# (For some simulations, it may make more sense
# to do this as part of the simulation driver.)
#------------------------------------------------------

calc_performance <- function(results) {
  
  performance_measures <- 
    results %>%
    filter(!is.na(p_val)) %>%
    group_by(cov_test, test) %>%
    summarize(K = n(),
              rej_rate_05 = mean(p_val < .05),
              mcse_05 = sqrt((rej_rate_05 * (1 - rej_rate_05))/K),
              rej_rate_01 = mean(p_val < .01),
              mcse_01 = sqrt((rej_rate_01 * (1 - rej_rate_01))/K),
              rej_rate_10 = mean(p_val < .10),
              mcse_10 = sqrt((rej_rate_10 * (1 - rej_rate_10))/K))

  return(performance_measures)
}
