library(tidyverse)
library(broom)
library(robumeta)
library(janitor)
library(xtable)
library(clubSandwich)
library(purrr)

load("../data/tsl_dat_20.RData")
source("2_estimation_study_1.R")

robu_comp_tsl <- robu(delta ~ g2age + dv, 
                      studynum = study, 
                      var.eff.size = v,
                      small = TRUE,
                      data = tsl_dat)

cr2_mat <- vcovCR(robu_comp_tsl, type = "CR2")

X <- model.matrix(as.formula("delta ~ g2age + dv"), data = tsl_dat)
y <- tsl_dat$delta
v <- tsl_dat$v
cluster <- tsl_dat$study


robu_comp_tsl

check_tsl <- robu_handmade(X = X, 
                           y = y, 
                           v = v, 
                           cluster = cluster, 
                           calc_vcov = "CR2")

check_tsl$vcov

near(check_tsl$vcov, cr2_mat)
near(check_tsl$coefficients, robu_comp_tsl$reg_table$b.r)


B_j <- attr(vcovCR(check_tsl, cluster = cluster, type = "CR2", inverse_var = TRUE), "adjustments")
B_j_club <- attr(vcovCR(robu_comp_tsl, type = "CR2"), "adjustments")

near(B_j[[15]], B_j_club[[15]])

naive_res <- Wald_test(check_tsl, 
                       constraints = constrain_zero(3:7),
                       vcov = cr2_mat,
                       test = "Naive-F", 
                       tidy = TRUE) %>%
  dplyr::select(`Naive-F` = p_val)




htz_res <- Wald_test(check_tsl, 
                     constraints = constrain_zero(3:7),
                     vcov = cr2_mat,
                     test = "HTZ", 
                     tidy = TRUE) %>%
  mutate(p_val = ifelse(Fstat < 0, 1, p_val)) %>%  # added this to fix the htz p val issue with negative F stat
  dplyr::select(HTZ = p_val)


# check cwb ---------------------------------------------------------------

tsl_dat$g <- tsl_dat$delta
tsl_dat$var_g <- tsl_dat$v

calculate_F <- function(beta, vcov, constraints, p = 7, test){
  
  C_mat <- diag(1L, nrow = p)[constraints,,drop = FALSE]    
  
  Q <- as.numeric(t(C_mat %*% beta) %*% chol2inv(chol(C_mat %*% vcov %*% t(C_mat))) %*% C_mat %*% beta)
  q <- nrow(C_mat)
  
  Fstat <- Q/q
  
  test_res <- tibble(test = test, 
                     Fstat = Fstat)
  
  return(test_res)
  
}

set.seed(7232020)

system.time(boot_tsl <- cwb(null_model = "g ~ g2age", 
    indices_test = 3:7, 
    R = 999, 
    full_form = "g2age + dv", 
    dat = tsl_dat))

save(boot_tsl, file = "../data/tsl_res_1108.RData")


Wald_test(robu_comp_tsl, 
          constraints = constrain_zero(3:7),
          test = "HTZ",
          vcov = "CR2")
