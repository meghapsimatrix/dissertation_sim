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



# Fit full model on data --------------------------------------------------
full_formula <- full_form
full_model <- robu(as.formula(paste("g ~ ", full_formula)), 
                   studynum = study, 
                   var.eff.size = var_g,
                   small = FALSE,
                   data = meta_data)

# get cov matrices --------------------------------------------------------
cov_mat_cr1 <- vcovCR(full_model, type = "CR1")
cov_mat_cr2 <- vcovCR(full_model, type = "CR2")

# get naive and htz -------------------------------------------------------
names(test_dat$indices_test) <- test_dat$cov_test
naive_res <- Wald_test(full_model, 
                       constraints = constrain_zero(test_dat$indices_test),
                       vcov = cov_mat_cr1,
                       test = "Naive-F", 
                       tidy = TRUE) %>%
  select(`Naive-F` = p_val)


htz_res <- Wald_test(full_model, 
                     constraints = constrain_zero(test_dat$indices_test),
                     vcov = cov_mat_cr2,
                     test = "HTZ", 
                     tidy = TRUE) %>%
  select(HTZ = p_val)


# cwb ---------------------------------------------------------------------
cwb_params <- test_dat %>%
  select(null_model, indices_test)

# if i don't put data and R and full_mod_form as default something goes wrong
# need to figure out how to do R 
# 490 seconds for 10 
# 2179.152  seconds for 31
system.time(boot_res <- pmap_dfr(cwb_params, cwb))
boot_res

save(boot_res, file = "../data/boot_res_1021.RData")

res <- 
  bind_cols(naive_res, htz_res, boot_res) %>%
  bind_cols(test_dat %>% select(cov_test, contrasts)) %>%
  gather(test, p_val, -c(cov_test, contrasts))

calc_performance(res)
