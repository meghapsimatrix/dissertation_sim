#------------------------------------------------------
# Model-fitting/estimation/testing functions
#------------------------------------------------------


# run the cwb -------------------------------------------------------------

change_to_mat <- function(res) {
  
  as.matrix(res)
  
}

mult_mat <- function(x, y) {
  
  as.vector(x %*% y)
  
}


extract_stats <- function(mod, C, vcov_mat, method) {
  
  Wald_test(mod, constraints = C, vcov = vcov_mat, test = "Naive-F") %>%
    as_tibble() %>%
    mutate(test = method) %>%
    select(test, Fstat)
}




cwb <- function(dat, null_model, R, boot_seed, indices_test) {
  
  null_mod <- robu(as.formula(null_model), 
                   studynum = study, 
                   var.eff.size = var_g,
                   small = FALSE,
                   data = dat)
  
  # residuals and transformed residuals -------------------------------------
  
  dat$res <- clubSandwich:::residuals_CS.robu(null_mod)
  dat$pred <- with(dat, g - res)
  split_res <- split(dat$res, dat$study)
  e_tilde_j <- map(split_res, change_to_mat)
  B_j <- attr(vcovCR(null_mod, type = "CR2"), "adjustments")
  dat$t_res <- unlist(pmap(list(B_j, e_tilde_j), mult_mat))
  
  
  # Rademacher weights ------------------------------------------------------
  
  num_studies <- unique(dat$study)
  k_j <- as.numeric(table(dat$study))
  
  set.seed(boot_seed)
  
  system.time(
    
    bootstraps <- rerun(.n = R, {
      
      wts <- sample(c(-1, 1), size = length(num_studies), replace = TRUE)
      dat$eta <- rep(wts, k_j)
      dat$new_t <- with(dat, pred + res * eta)
      dat$new_t_adj <- with(dat, pred + t_res * eta)
      
      full_mod <- fit_mod("new_t ~ X1 + X2 + X3 + X4 + X5", dat)
      full_mod_adj <- fit_mod("new_t_adj ~ X1 + X2 + X3 + X4 + X5", dat)
      
      cov_mat <- vcovCR(full_mod, type = "CR1")
      cov_mat_adj <- vcovCR(full_mod_adj, type = "CR1")
      
      res <- extract_stats(full_mod, constrain_zero(indices_test), cov_mat, "CWB")
      res_adj <- extract_stats(full_mod_adj, constrain_zero(indices_test), cov_mat_adj, "CWB Adjusted")
      
      bind_rows(res, res_adj)
      
    }) %>%
      bind_rows()
    
  )
  
  org_F <- Wald_test(full_mod, 
                     constraints = constrain_zero(indices_test), 
                     vcov = vcovCR(full_mod, type = "CR1"),
                     test = "Naive-F") %>%
    pull(Fstat)
  
  
  p_boot <- 
    bootstraps %>%
    group_by(test) %>%
    summarize(p_val = mean(Fstat > org_F)) %>%
    ungroup() %>%
    spread(test, p_val)
  
  
  return(p_boot)
  
}

