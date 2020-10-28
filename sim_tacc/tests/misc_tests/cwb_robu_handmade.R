library(dplyr)
library(purrr)
library(mvtnorm)
library(robumeta)
library(clubSandwich)
library(tidyr)


load("../data/meta_data_practice.RData")
load("data/to_test.RData")
load("data/design_mat.Rdata")

glimpse(to_test)

source("2_estimation_study_1.R")
source("3_performance_criteria.R")
source("1_data_gen_study_1.R")


test_dat <- to_test
rm(to_test)
full_form <- "X1 + X2 + X3 + X4 + X5"
R <- 399


dat <- meta_data
null_model <- test_dat$null_model[[1]]

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
  

    
    wts <- sample(c(-1, 1), size = length(num_studies), replace = TRUE)
    dat$eta <- rep(wts, k_j)
    new_t <- with(dat, pred + res * eta)
    new_t_adj <- with(dat, pred + t_res * eta)
    
    y <- dat$g
    v <- dat$var_g
    X <- model.matrix(g ~ X1 + X2 + X3 + X4 + X5, data = dat)
    cluster <- dat$study
    
    handmade_fit <- robu_handmade(X = X, y = y, v = v, cluster = cluster)
    full_mod <- update_robu(handmade_fit, y = new_t)
    full_mod_adj <- update_robu(handmade_fit, y = new_t_adj)
    
    cov_mat <- vcovCR(full_mod, type = "CR1")
    cov_mat_adj <- vcovCR(full_mod_adj, type = "CR1")
    
    res <- calculate_F(beta = full_mod$reg_table$b.r, vcov = cov_mat, constraints = indices_test, test = "CWB")
    res_adj <- calculate_F(beta = full_mod_adj$reg_table$b.r, vcov = cov_mat_adj, constraints = indices_test, test = "CWB Adjusted")
    
    bind_rows(res, res_adj)
    

  
  
  
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


