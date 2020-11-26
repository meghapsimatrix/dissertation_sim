library(dplyr)
library(purrr)
library(mvtnorm)
library(clubSandwich)
library(tidyr)


load("../data/meta_data_practice_2.Rdata")


source("2_estimation_study_2.R")
source("3_performance_criteria.R")
source("1_data_gen_study_2.R")



full_form <- "X1_B + X1_C + X1_D + X1_E"
R <- 399


#set.seed(11032020)

# meta_data <-
#   generate_rmeta(m = 10,
#                  tau = 0.1,
#                  rho = 0.8,
#                  covs = design_mat,
#                  beta_type = "A")


# Fit full model on data --------------------------------------------------
y <- meta_data$g
v <- meta_data$var_g
X <- model.matrix(as.formula(paste("g ~ ", full_form)), data = meta_data)
cluster <- meta_data$study

full_model <- robu_handmade(X = X, y = y, v = v, cluster = cluster)

full_model

# get cov matrices --------------------------------------------------------
cov_mat_cr1 <- vcovCR(full_model, type = "CR1", cluster = cluster)
cov_mat_cr2 <- vcovCR(full_model, type = "CR2", cluster = cluster)

# get naive and htz -------------------------------------------------------
naive_res <- Wald_test(full_model, 
                       constraints = constrain_zero(2:5),
                       vcov = cov_mat_cr1,
                       test = "Naive-F", 
                       tidy = TRUE) %>%
  select(`Naive-F` = p_val)


htz_res <- Wald_test(full_model, 
                     constraints = constrain_zero(2:5),
                     vcov = cov_mat_cr2,
                     test = "HTZ", 
                     tidy = TRUE) %>%
  select(HTZ = p_val)

# the last one is NA?
htz_res %>% View()




# cwb ---------------------------------------------------------------------

full_model$coefficients

system.time(res <- cwb(null_model = "g ~ 1", 
                indices_test = 2:5, 
                full_form = full_form, 
                R = R, dat = meta_data))
