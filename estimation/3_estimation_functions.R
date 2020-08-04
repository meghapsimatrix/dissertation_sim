

estimate <- function(full_formula, null_formula, dat){
  
  full_mod <- robu(as.formula(full_formula), 
                   studynum = study, 
                   var.eff.size = v,
                   small = FALSE,
                   data = dat)
  
  # Single coefficient ------------------------------------------------------
  
  cov_mat_cr0 <- vcovCR(full_mod, type = "CR0")
  cov_mat_cr2 <- vcovCR(full_mod, type = "CR2")
  
  res_s_naive_b2 <- Wald_test(full_mod, 
                              constraints = constrain_zero(2), 
                              vcov = cov_mat_cr0,
                              test = "Naive-F")
  
  res_s_htz_b2  <- Wald_test(full_mod, 
                             constraints = constrain_zero(2), 
                             vcov = cov_mat_cr2,
                             test = "HTZ")
  
  # multiple contrast hypothesis --------------------------------------------
  
  res_mch_naive <- Wald_test(full_mod, 
                             constraints = constrain_zero(2:6), 
                             vcov = cov_mat_cr0, 
                             test = "Naive-F")
  
  res_mch_htz <- Wald_test(full_mod, 
                           constraints = constrain_zero(2:6), 
                           vcov = cov_mat_cr2,
                           test = "HTZ")
  
}


estimate_cwb <- function(dat){
  
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
  
  res_single <- extract_stats(full_mod, constrain_zero(2), cov_mat, "CWB", "age")
  res_single_adj <- extract_stats(full_mod_adj, constrain_zero(2), cov_mat_adj, "CWB Adjusted", "age")
  res_mch <- extract_stats(full_mod, constrain_zero(3:7), cov_mat, "CWB", "mch")
  res_mch_adj <- extract_stats(full_mod_adj, constrain_zero(3:7), cov_mat_adj, "CWB Adjusted", "mch")
  
  res <- bind_rows(res_single, res_single_adj, res_mch, res_mch_adj)
  
  return(res)
  
}