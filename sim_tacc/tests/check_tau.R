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

cov_mat_cr1 <- vcovCR(full_model, type = "CR1", cluster = cluster)
cov_mat_cr2 <- vcovCR(full_model, type = "CR2", cluster = cluster)


htz_res <- Wald_test(full_model, 
                     constraints = constrain_zero(test_dat$indices_test),
                     vcov = cov_mat_cr2,
                     test = "HTZ", 
                     tidy = TRUE) 


%>%
  mutate(p_val = ifelse(Fstat < 0, 1, p_val)) %>%  # added this to fix the htz p val issue with negative F stat
  dplyr::select(HTZ = p_val)




k_j <- as.numeric(table(cluster))
sigma_sq_j <- tapply(v, cluster, mean)

m <- length(k_j)

w_tilde <- 1 / (k_j * sigma_sq_j)
w_tilde
w_tilde_j <- rep(w_tilde, k_j)

mod_prelim <- lm.wfit(x = X, y = y, w = w_tilde_j)

res <- residuals(mod_prelim) 

res

# calculate weighted residual sum of squares
QE <- sum(w_tilde_j * res^2)
QE

# split the design matrix by study
X_j <- mat_split(X, cluster)
X_j

# Create M tilde ----------------------------------------------------------
M_tilde <- chol2inv(chol(crossprod(X, w_tilde_j * X)))
p_j <- lapply(X_j, colSums)

# trace products ----------------------------------------------------------

# the first B for the trace product numerator
w_over_k <- w_tilde / k_j
B_num_all_1 <- Map(function(w, X) w * crossprod(X), w = w_over_k, X = X_j)
B_num_1 <- reduce(B_num_all_1, `+`)

# the second B for the trace product numerator
XJX_j <- lapply(p_j, tcrossprod)
B_num_all_2 <- Map(function(w, xjx) w * xjx, w = w_over_k, xjx = XJX_j)
B_num_2 <- reduce(B_num_all_2, `+`)

# the B for the trade product denominator
B_den_all <- Map(function(w, xjx) w^2 * xjx, w = w_tilde, xjx = XJX_j)
B_den <- reduce(B_den_all, `+`)

num_minus <- m - (1 - rho) * trace_product(M_tilde, B_num_1) - rho * trace_product(M_tilde, B_num_2)
den <- sum(k_j * w_tilde) - trace_product(M_tilde, B_den)

tau_sq <- (QE - num_minus) / den
tau_sq

# CE weights
w_j <- 1 / (k_j * (sigma_sq_j + tau_sq))
w_j
w_ij <- rep(w_j, k_j)

# fit WLS regression
mod_CE <- lm.wfit(x = X, y = y, w = w_ij)

res <- mod_CE[c("coefficients","residuals","fitted.values","weights")]
res$tau_sq <- tau_sq

res$X <- X
res$w_ij <- w_ij
res$cluster <- cluster
res$k_j <- k_j
res$sigma_sq_j <- sigma_sq_j
res$w_tilde_j <- w_tilde_j
res$num_minus <- num_minus
res$den <- den
res$nobs <- sum(k_j)

