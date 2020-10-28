library(simhelpers)


calc_performance <- function(sim_dat){
  
  rr_dat <- sim_dat %>%
    group_by(m, tau, rho, beta, test, indices) %>%
    group_modify(~ calc_rejection(.x, p_values = p_val))
  
  return(rr_dat)
  
}