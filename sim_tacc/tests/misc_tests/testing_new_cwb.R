

library(tidyverse)
library(broom)
library(robumeta)
library(janitor)
library(xtable)
library(clubSandwich)

load("../data/tsl_dat_20.RData")

full_model <- robu(delta ~ g2age + dv, 
                      studynum = study, 
                      var.eff.size = v,
                      small = TRUE,
                      data = tsl_dat)

full_model

cov_mat_first <- vcovCR(full_model, type = "CR1")

# run the cwb -------------------------------------------------------------

change_to_mat <- function(res) {
  
  as.matrix(res)
  
}

mult_mat <- function(x, y) {
  
  as.vector(x %*% y)
  
}

calculate_F <- function(beta, vcov, constraints, p = 7, test){
  
  C_mat <- diag(1L, nrow = p)[constraints,,drop = FALSE]    
  
  Q <- as.numeric(t(C_mat %*% beta) %*% chol2inv(chol(C_mat %*% vcov %*% t(C_mat))) %*% C_mat %*% beta)
  q <- nrow(C_mat)
  
  Fstat <- Q/q
  
  test_res <- tibble(test = test, 
                     Fstat = Fstat)
  
  return(test_res)
  
}


cwb <- function(null_model, 
                indices_test, 
                R, 
                full_form, 
                full_mod_org, 
                cov_mat_org, 
                dat) {
  
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
  
  bootstraps <- rerun(.n = R, {
    
    wts <- sample(c(-1, 1), size = length(num_studies), replace = TRUE)
    dat$eta <- rep(wts, k_j)
    dat$new_t <- with(dat, pred + res * eta)
    dat$new_t_adj <- with(dat, pred + t_res * eta)
    
    full_mod <- robu(as.formula(paste("new_t ~", full_form)), 
                     studynum = study, 
                     var.eff.size = var_g,
                     small = FALSE,
                     data = dat)
    
    full_mod_adj <- robu(as.formula(paste("new_t_adj ~", full_form)), 
                         studynum = study, 
                         var.eff.size = var_g,
                         small = FALSE,
                         data = dat)
    
    
    cov_mat <- vcovCR(full_mod, type = "CR1")
    cov_mat_adj <- vcovCR(full_mod_adj, type = "CR1")
    
    res <- calculate_F(beta = full_mod$reg_table$b.r, vcov = cov_mat, constraints = indices_test, test = "CWB")
    res_adj <- calculate_F(beta = full_mod_adj$reg_table$b.r, vcov = cov_mat_adj, constraints = indices_test, test = "CWB Adjusted")
    
    bind_rows(res, res_adj)
    
  }) %>%
    bind_rows()
  
  
  
  org_F <- calculate_F(beta = full_mod_org$reg_table$b.r, 
                       vcov = cov_mat_org, 
                       constraints = indices_test, 
                       test = "Naive") %>%
    pull(Fstat)
  
  
  p_boot <- 
    bootstraps %>%
    group_by(test) %>%
    summarize(p_val = mean(Fstat > org_F)) %>%
    ungroup() %>%
    spread(test, p_val)
  
  
  return(p_boot)
  
}


set.seed(7232020)

tsl_dat <- tsl_dat %>%
  mutate(var_g = v, 
         g = delta)

system.time(res <- cwb(null_model = "g ~ g2age", indices_test = 3:7, R = 999, 
    full_form = "g2age + dv",
    full_mod_org = full_model,
    cov_mat_org = cov_mat_first,
    dat = tsl_dat))





res_mch_htz <- Wald_test(full_model, 
                         constraints = constrain_zero(3:7), 
                         vcov = "CR2",
                         test = "HTZ")

res_mch_htz
