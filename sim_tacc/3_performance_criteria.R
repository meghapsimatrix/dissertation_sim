#------------------------------------------------------
# Calculate performance measures
# (For some simulations, it may make more sense
# to do this as part of the simulation driver.)
#------------------------------------------------------

calc_performance <- function(results, alpha) {
  
  performance_measures  <- results %>%
    filter(!is.na(p_val)) %>%
    group_by(test) %>%
    summarize(K = n(),
              rej_rate = mean(p_val < alpha),
              mcse = sqrt((rej_rate * (1 - rej_rate))/K))
  
  return(performance_measures)
}
