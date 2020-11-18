extract_stats <- function(mod, constraints, vcov_mat, method, var){
  
  Wald_test(mod, constraints = constraints, vcov = vcov_mat, test = "Naive-F") %>%
    as_tibble() %>%
    mutate(type = method,
           constraint = var)
}

cwb <- function(dat, single){
  
  num_studies <- unique(dat$study)
  wts <- sample(c(-1, 1), size = length(num_studies), replace = TRUE)
  k_j <- as.numeric(table(dat$study))
  
  dat$eta <- rep(wts, k_j)
  
  dat$new_t <- with(dat, pred_null + res_null * eta)
  dat$new_t_adj <- with(dat, pred_null + t_res * eta)
  
  
  full_mod <- robu(new_t ~ g2age + dv, 
                   studynum = study, 
                   var.eff.size = v,
                   small = FALSE,
                   data = dat)
  
  full_mod_adj <- robu(new_t_adj ~ g2age + dv, 
                       studynum = study, 
                       var.eff.size = v,
                       small = FALSE,
                       data = dat)
  
  cov_mat <- vcovCR(full_mod, type = "CR2")
  cov_mat_adj <- vcovCR(full_mod_adj, type = "CR2")
  
  if(single == TRUE){
    res <- extract_stats(full_mod, constrain_zero(2), cov_mat, "CWB", "age")
    res_adj <- extract_stats(full_mod_adj, constrain_zero(2), cov_mat_adj, "CWB Adjusted", "age")
  }
  else{
    res <- extract_stats(full_mod, constrain_zero(3:7), cov_mat, "CWB", "mch")
    res_adj <- extract_stats(full_mod_adj, constrain_zero(3:7), cov_mat_adj, "CWB Adjusted", "mch")
  }
  
  res <- bind_rows(res, res_adj)
  
  return(res)
  
}