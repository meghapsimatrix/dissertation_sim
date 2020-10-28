library(dplyr)
library(purrr)
library(mvtnorm)
library(robumeta)
library(clubSandwich)
library(tidyr)
library(stringr)

# Tipton Pusto design matrix cleaned - clean_design_mat.R
load("sim_tacc/data/design_mat.Rdata")
load("sim_tacc/data/to_test.Rdata")

#-----------------------------------------------------------
# Source the functions
#-----------------------------------------------------------

source("sim_tacc/1_data_gen_study_1.R")
source("sim_tacc/2_estimation_study_1.R")
source("sim_tacc/3_performance_criteria.R")
source("homemade/handmade_robu.R")

meta_data <- 
  generate_rmeta(m = 80, 
                 tau = 0.4, 
                 rho = 0.8, 
                 covs = design_mat,
                 beta_type = "A")

robu_fit <- robu(g ~ X1 + X2 + X3 + X4 + X5, 
                 data = meta_data, studynum = study, 
                 var.eff.size = var_g)

y <- meta_data$g
v <- meta_data$var_g
X <- model.matrix(g ~ X1 + X2 + X3 + X4 + X5, data = meta_data)
cluster <- meta_data$study

handmade_fit <- robu_handmade(X = X, y = y, v = v, cluster = cluster)
update_fit <- update_robu(handmade_fit, y = y)


full_mod <-  robu_handmade(X = X, y = y, v = v, cluster = cluster, calc_vcov = "CR2")
full_mod$vcov


B_j <- attr(vcovCR(handmade_fit, cluster = cluster, type = "CR2", inverse_var = TRUE), "adjustments")
B_j_club <- attr(vcovCR(robu_fit, type = "CR2"), "adjustments")

near(B_j[[60]], B_j_club[[60]])

handmade_fit$cluster
update_fit$residuals


update_fit$vcov
from_club <- vcovCR(robu_fit, type = "CR0")
near(update_fit$vcov, from_club)

# vcovCR now works with handmade_robu() models
CR0 <- vcovCR(handmade_fit, cluster = cluster, type = "CR0")
CR0
update_fit$vcov
near(update_fit$vcov, CR0)

# be sure to set inverse_var = TRUE when calculating CR2
CR2 <- vcovCR(handmade_fit, cluster = cluster, type = "CR2", inverse_var = TRUE)
CR2
robu_fit$VR.r
from_club_CR2 <- vcovCR(robu_fit, type = "CR2")
near(CR2, robu_fit$VR.r)
near(CR2, from_club_CR2)
near(full_mod$vcov, from_club_CR2)


# tau-squared
c(handmade_fit$tau_sq, update_fit$tau_sq, as.numeric(robu_fit$mod_info$tau.sq))
all.equal(handmade_fit$tau_sq, as.numeric(robu_fit$mod_info$tau.sq))
all.equal(handmade_fit$tau_sq, update_fit$tau_sq)

# beta
cbind(handmade_fit$coefficients,
      update_fit$coefficients,
      robu_fit$reg_table$b.r)
all.equal(handmade_fit$coefficients, robu_fit$reg_table$b.r, check.attributes = FALSE)
all.equal(handmade_fit$coefficients, update_fit$coefficients)

library(microbenchmark)

microbenchmark(
  robu = robu(g ~ X1 + X2 + X3 + X4 + X5, 
              data = meta_data, studynum = study, 
              var.eff.size = var_g, small = FALSE),
  hand = robu_handmade(X = X, y = y, v = v, cluster = cluster),
  up = update_robu(handmade_fit, y = y)
)

microbenchmark(
  robu = robu(g ~ X1 + X2 + X3 + X4 + X5, 
              data = meta_data, studynum = study, 
              var.eff.size = var_g, small = TRUE),
  hand = robu_handmade(X = X, y = y, v = v, cluster = cluster, calc_vcov = "CR2")
)
