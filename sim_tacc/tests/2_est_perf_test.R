library(dplyr)
library(purrr)
library(mvtnorm)
library(robumeta)
library(clubSandwich)
library(tidyr)


load("../data/meta_data_practice.RData")
load("data/to_test.RData")

glimpse(to_test)



source("2_estimation_study_1.R")
source("3_performance_criteria.R")


test_dat <- to_test
rm(to_test)
full_form <- "g ~ X1 + X2 + X3 + X4 + X5"

# Fit full model on data --------------------------------------------------
full_formula <- full_form
full_model <- robu(as.formula(full_formula), 
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

naive_res <- map_dfr(test_dat$indices_test,
                     estimate_wald, 
                     model = full_model, 
                     cov_mat = cov_mat_cr1, 
                     test = "Naive-F")


htz_res <-  Wald_test(full_model, 
                      constraints = constrain_zero(test_dat$indices_test),
                      vcov = cov_mat_cr2,
                      test = "HTZ", 
                      tidy = TRUE) %>%
  select(HTZ = p_val)


# cwb ---------------------------------------------------------------------
set.seed(10152020)

test_dat <- test_dat %>%
  mutate(boot_seed = round(runif(1) * 2^30) + 1:n())

cwb_params <- test_dat %>%
  select(null_model, R, boot_seed, indices_test)

boot_res <- pmap_dfr(cwb_params, cwb)

res <- 
  bind_cols(naive_res, htz_res, boot_res) %>%
  bind_cols(test_dat %>% select(cov_test, contrasts)) %>%
  gather(test, p_val, -c(cov_test, contrasts))





