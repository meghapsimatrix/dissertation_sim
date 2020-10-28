library(robumeta)
library(tidyverse)
library(clubSandwich)

load("data/meta_data_practice.RData")
load("sim_tacc/data/to_test.RData")

dat <- meta_data


full_formula <- "g ~ X1 + X2 + X3 + X4 + X5"
full_model <- robu(as.formula(full_formula), 
                   studynum = study, 
                   var.eff.size = var_g,
                   small = FALSE,
                   data = dat)

estimate_wald <- function(model, indices_test, cov_mat, test){
  
  res <- Wald_test(model, 
                   constraints = constrain_zero(indices_test), 
                   vcov = cov_mat,
                   test = test) %>%
    select(test, p_val)
  
  return(res)
}



# run wald ----------------------------------------------------------------

cov_mat_cr1 <- vcovCR(full_model, type = "CR1")
cov_mat_cr2 <- vcovCR(full_model, type = "CR2")

system.time(naive_res <- to_test %>%
  mutate(res = map(indices_test, estimate_wald, 
                       model = full_model, 
                       cov_mat = cov_mat_cr1, 
                       test = "Naive-F")) %>%
  unnest(cols = res))


system.time(htz_res <- to_test %>%
  mutate(res = map(indices_test, 
                   estimate_wald, 
                   model = full_model, 
                   cov_mat = cov_mat_cr2, 
                   test = "HTZ")) %>%
  unnest(cols = res))





