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


# Fit full model on data --------------------------------------------------
full_model <- robu(as.formula(paste("g ~ ", full_form)), 
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
  select(null_model, indices_test) %>%
  mutate(R = R,
         full_form = full_form)

system.time(boot_res <- pmap_dfr(cwb_params[1:2, ], 
                                 .f = cwb, 
                                 dat = meta_data))
small_boot <- boot_res
small_boot


set.seed(10282020)

system.time(boot_res <- pmap_dfr(cwb_params, 
                                 .f = cwb, 
                                 dat = meta_data))
boot_res

# 1502.069  on 1026 afternoon
# 140.722   on 1027 night
save(boot_res, file = "../data/boot_res_1028.RData")

res <- 
  bind_cols(naive_res, htz_res, boot_res) %>%
  bind_cols(test_dat %>% select(cov_test, contrasts)) %>%
  gather(test, p_val, -c(cov_test, contrasts))

calc_performance(res)
