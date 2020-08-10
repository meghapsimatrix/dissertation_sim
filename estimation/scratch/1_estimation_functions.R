

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


