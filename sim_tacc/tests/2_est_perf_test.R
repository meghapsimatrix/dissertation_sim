library(dplyr)
library(purrr)
library(mvtnorm)
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

# the last one is NA?
htz_res %>% View()

Wald_test(full_model, 
          constraints = constrain_zero(c(2, 3, 4, 5, 6)),
          vcov = cov_mat_cr2,
          test = "HTZ")

# Error in symnum(x$p_val, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001,  : 
# must have 2 'symbols' for logical 'x' argument


# cwb ---------------------------------------------------------------------
cwb_params <- test_dat %>%
  select(null_model, indices_test)

system.time(boot_res <- pmap_dfr(cwb_params[1:2, ], 
                                 .f = cwb, 
                                 dat = meta_data,
                                 R = R,
                                 full_form = full_form))
small_boot <- boot_res
small_boot


set.seed(10282020)

system.time(boot_res <- pmap_dfr(cwb_params, 
                                 .f = cwb, 
                                 dat = meta_data,
                                 R = R,
                                 full_form = full_form))
boot_res

# 1502.069  on 1026 afternoon
# 140.722   on 1027 night
# 171.116  on 1030
save(boot_res, file = "../data/boot_res_1108.RData")

boot_res_today <- boot_res

load("../data/boot_res_1030.RData")

near(boot_res_today, boot_res)

res <- 
  bind_cols(naive_res, htz_res, boot_res) %>%
  bind_cols(test_dat %>% select(cov_test, contrasts)) %>%
  gather(test, p_val, -c(cov_test, contrasts))

results <- calc_performance(res)



# handmade tests ----------------------------------------------------------
library(robumeta)

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

