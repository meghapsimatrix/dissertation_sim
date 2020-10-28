library(dplyr)
library(purrr)
library(mvtnorm)
library(robumeta)
library(clubSandwich)
library(tidyr)
library(stringr)

# Tipton Pusto design matrix cleaned - clean_design_mat.R
load("data/design_mat.Rdata")
load("data/to_test.Rdata")

#-----------------------------------------------------------
# Source the functions
#-----------------------------------------------------------

source("1_data_gen_study_1.R")
source("2_estimation_study_1.R")
source("3_performance_criteria.R")

m <- 10
tau <- .1
rho <- .5
design_matrix <- design_mat
beta_type <- "A"
full_form <- "X1 + X2 + X3 + X4 + X5"
test_dat <- to_test
R <- 20
iterations <- 2


results <-
  rerun(iterations, {
    
    # generate data ------------------------------------------------------------
    meta_data <- generate_rmeta(m = m, 
                                tau = tau, 
                                rho = rho,
                                covs = design_matrix,
                                beta_type = beta_type)
    
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
    
    # if i don't put data and R and full_mod_form as default something goes wrong
    # need to figure out how to do R 
    boot_res <- pmap_dfr(cwb_params, cwb, dat = meta_data)
    
    res <- 
      bind_cols(naive_res, htz_res, boot_res) %>%
      bind_cols(test_dat %>% dplyr::select(cov_test)) %>%
      gather(test, p_val, -c(cov_test))
    
  })  %>%
  bind_rows()

calc_performance(results)

