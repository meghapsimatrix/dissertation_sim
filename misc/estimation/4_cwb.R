library(robumeta)
library(tidyverse)
library(clubSandwich)

load("data/meta_data_practice.RData")
load("sim_tacc/data/to_test.RData")


full_formula <- "g ~ X1 + X2 + X3 + X4 + X5"
full_model <- robu(as.formula(full_formula), 
                   studynum = study, 
                   var.eff.size = var_g,
                   small = FALSE,
                   data = meta_data)


# cwb ---------------------------------------------------------------------

fit_mod <- function(equation, dat = meta_data) {
  
  robu(as.formula(equation), 
       studynum = study, 
       var.eff.size = var_g,
       small = FALSE,
       data = dat)
  
}


system.time(null_mods <- map(to_test$null_model, fit_mod))




# run the cwb -------------------------------------------------------------

change_to_mat <- function(res){
  
  as.matrix(res)
  
}

mult_mat <- function(x, y) {
  
  as.vector(x %*% y)
  
}


extract_stats <- function(mod, C, vcov_mat, method){
  
  Wald_test(mod, constraints = C, vcov = vcov_mat, test = "Naive-F") %>%
    as_tibble() %>%
    mutate(test = method) %>%
    select(test, Fstat)
}



cwb <- function(dat = meta_data, null_mod, full_mod = full_model, indices_test) {
  
  
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
  
  set.seed(8062020)
  
  system.time(
    
    bootstraps <- rerun(.n = 399, {
      
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
  
  
  p_boot <- bootstraps %>%
    group_by(test) %>%
    summarize(p_val = mean(Fstat > org_F)) %>%
    ungroup() 
  
  
  return(p_boot)
  
}


# Similar to HTZ and Naive F

system.time(res <- cwb(dat = meta_data, 
                       null_mod = to_test$null_mod[[25]], 
                       indices_test = to_test$indices_test[[25]]))

Wald_test(full_model, constraints = constrain_zero(c(4, 5, 6)), 
          vcov = vcovCR(full_model, type = "CR2"), test = "HTZ")

Wald_test(full_model, constraints = constrain_zero(c(4, 5, 6)), 
          vcov = vcovCR(full_model, type = "CR1"), test = "Naive-F")

res

map2(null_mods[1:5], indices_test[1:5], cwb)


cwb_params <- to_test %>%
  mutate(null_mod = null_mods) %>%
  select(null_mod, indices_test)

system.time(cwb_res <- pmap_dfr(cwb_params[1:2, ], cwb))


