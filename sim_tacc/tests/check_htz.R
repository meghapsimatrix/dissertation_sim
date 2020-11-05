library(dplyr)
library(purrr)
library(mvtnorm)
library(robumeta)
library(clubSandwich)
library(tidyr)


load("data/design_mat.Rdata")
load("data/to_test.RData")

#-----------------------------------------------------------
# Source the functions
#-----------------------------------------------------------

source("1_data_gen_study_1.R")
source("2_estimation_study_1.R")

set.seed(509985477)


m <- 10
tau <- 0.1 
rho <- 0.8 
design_matrix <- design_mat
beta_type <- "A"
test_dat <- to_test
rm(to_test)
full_form <- "X1 + X2 + X3 + X4 + X5"
R <- 50


results <-
  rerun(2, {
    
    # generate data ------------------------------------------------------------
    meta_data <- generate_rmeta(m = m, 
                                tau = tau, 
                                rho = rho,
                                covs = design_matrix,
                                beta_type = beta_type)
    
    
  }) 






meta_data <- results[[2]]
library(robumeta)

robu_fit <- robu(g ~ X1 + X2 + X3 + X4 + X5, 
                 data = meta_data, studynum = study, 
                 var.eff.size = var_g)

robu_fit$mod_info$tau.sq

# Fit full model on data --------------------------------------------------
y <- meta_data$g
v <- meta_data$var_g
X <- model.matrix(as.formula(paste("g ~ ", full_form)), data = meta_data)
cluster <- meta_data$study

full_model <- robu_handmade(X = X, y = y, v = v, cluster = cluster)

full_model$tau_sq


# get cov matrices --------------------------------------------------------
cov_mat_cr1 <- vcovCR(full_model, type = "CR1", cluster = cluster)
cov_mat_cr2 <- vcovCR(full_model, type = "CR2", cluster = cluster)

# get naive and htz -------------------------------------------------------
names(test_dat$indices_test) <- test_dat$cov_test
naive_res <- Wald_test(full_model, 
                       constraints = constrain_zero(test_dat$indices_test),
                       vcov = cov_mat_cr1,
                       test = "Naive-F", 
                       tidy = TRUE) %>%
  dplyr::select(`Naive-F` = p_val)


htz_res <- Wald_test(full_model, 
                     constraints = constrain_zero(test_dat$indices_test),
                     vcov = cov_mat_cr2,
                     test = "HTZ", 
                     tidy = TRUE) %>%
  dplyr::select(HTZ = p_val)

htz_res %>%
  View()


htz_res <- Wald_test(full_model, 
                     constraints = constrain_zero(test_dat$indices_test),
                     vcov = cov_mat_cr2,
                     test = "HTZ", 
                     tidy = TRUE) %>%
  mutate(p_val = ifelse(Fstat < 0, 1, p_val)) %>%  # added this to fix the htz p val issue with negative F stat
  dplyr::select(HTZ = p_val)

htz_res %>%
  View()
