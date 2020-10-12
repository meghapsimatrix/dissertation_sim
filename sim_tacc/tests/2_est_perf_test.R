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


# Fit full model on data --------------------------------------------------
full_formula <- "g ~ X1 + X2 + X3 + X4 + X5"
full_model <- robu(as.formula(full_formula), 
                   studynum = study, 
                   var.eff.size = var_g,
                   small = FALSE,
                   data = meta_data)


# get cov matrices --------------------------------------------------------

cov_mat_cr1 <- vcovCR(full_model, type = "CR1")
cov_mat_cr2 <- vcovCR(full_model, type = "CR2")

# get naive and htz -------------------------------------------------------

naive_res <- map_dfr(test_dat$indices_test,
                     estimate_wald, 
                     model = full_model, 
                     cov_mat = cov_mat_cr1, 
                     test = "Naive-F")

naive_res$`Naive-F`

# JEP: Another way to do the same thing

names(test_dat$indices_test) <- test_dat$cov_test
Wald_test(full_model, 
          constraints = constrain_zero(test_dat$indices_test),
          vcov = cov_mat_cr1,
          test = "Naive-F", 
          tidy = TRUE)


htz_res <-  map_dfr(test_dat$indices_test,
                    estimate_wald, 
                    model = full_model, 
                    cov_mat = cov_mat_cr2, 
                    test = "HTZ")
htz_res$HTZ

# JEP: Another way to do the same thing

Wald_test(full_model, 
          constraints = constrain_zero(test_dat$indices_test),
          vcov = cov_mat_cr2,
          test = "HTZ", 
          tidy = TRUE)


# cwb ---------------------------------------------------------------------

null_mods <- map(test_dat$null_model, fit_mod)

cwb_params <- test_dat %>%
  mutate(null_mod = null_mods) %>%
  select(null_mod, indices_test)

boot_res <- pmap_dfr(cwb_params[1:2, ], cwb)  # just doing two boots

boot_res_2 <- pmap_dfr(cwb_params[20, ], cwb) # checking mch

res <- bind_cols(naive_res, htz_res) %>%  # the one on run_sim includes boot_res, here just to check
  bind_cols(test_dat %>% select(cov_test, contrasts)) %>%
  gather(test, p_val, -c(cov_test, contrasts))


calc_performance(res)



# JEP scratch ---------------------------------------------------------
test <- 20
dat <- meta_data
null_mod <- null_mods[[test]]
full_mod <- full_model
indices_test <- test_dat$indices_test[[test]]


iterations <- 4
m <- 20
tau <- 0.1
rho <- 0.6
beta_type <- "A"
test_dat <- to_test[1:3,]
seed <- NULL
