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

y <- dat$g
v <- dat$var_g
X_null <- model.matrix(as.formula(null_model), data = dat)
X_full <- model.matrix(as.formula(paste("g ~ ", full_form)), data = dat)

cluster <- dat$study

null_mod <- robu_handmade(X = X_null, y = y, v = v, cluster = cluster)
full_mod_org <- robu_handmade(X = X_full, y = y, v = v, cluster = cluster, calc_vcov = "CR0")
cov_mat_org <- full_mod_org$vcov

# residuals and transformed residuals -------------------------------------

dat$res <- null_mod$residuals
dat$pred <- null_mod$fitted.values
split_res <- split(dat$res, dat$study)
e_tilde_j <- map(split_res, change_to_mat)
B_j <- attr(vcovCR(null_mod, cluster = cluster, type = "CR2", inverse_var = TRUE), "adjustments")
dat$t_res <- unlist(pmap(list(B_j, e_tilde_j), mult_mat))


# Rademacher weights ------------------------------------------------------

num_studies <- unique(dat$study)
k_j <- as.numeric(table(dat$study))

bootstraps <- rerun(.n = R, {
  
  wts <- sample(c(-1, 1), size = length(num_studies), replace = TRUE)
  dat$eta <- rep(wts, k_j)
  new_t <- with(dat, pred + res * eta)
  new_t_adj <- with(dat, pred + t_res * eta)
  
  
  
  full_mod_cwb <- update_robu(full_mod_org, y = new_t)
  full_mod_cwb_adj <- update_robu(full_mod_org, y = new_t_adj)
  
  cov_mat_cwb <- full_mod_cwb$vcov
  cov_mat_cwb_adj <- full_mod_cwb_adj$vcov
  
  res <- calculate_F(beta = full_mod_cwb$coefficients, vcov = cov_mat_cwb, constraints = indices_test, test = "CWB")
  res_adj <- calculate_F(beta = full_mod_cwb_adj$coefficients, vcov = cov_mat_cwb_adj, constraints = indices_test, test = "CWB Adjusted")
  
  bind_rows(res, res_adj)
  
}) %>%
  bind_rows()



org_F <- calculate_F(beta = full_mod_org$coefficients, 
                     vcov = cov_mat_org, 
                     constraints = indices_test, 
                     test = "Naive") %>%
  pull(Fstat)

