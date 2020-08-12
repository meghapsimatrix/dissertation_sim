library(tidyverse)
library(mvtnorm)

# my notation is different from the one in James's blog post btw :D 


# Generate one study  -----------------------------------------------------

# m is number of study
# k is number of effect sizes per study
# j is the jth study out of m studies
# N is number of participants per study
# Psi is the correlation (or correlation matrix) of the participant-level outcomes (conditional on true effect size delta). 

generate_rsmd <- function(delta, k, N, Psi) {
  
  # make sure delta is a vector
  delta_vec <- rep(delta, length.out = k)
  
  # create Psi matrix assuming equicorrelation
  if (!is.matrix(Psi)) Psi <- Psi + diag(1 - Psi, nrow = k) # cor matrix for 1 study
  
  # generate numerator of SMD 
  meandiff <- rmvnorm(n = 1, mean = delta_vec, sigma = (4/N) * Psi) 
  
  # covariance 
  cov_mat <- as.matrix(rWishart(n = 1, df = N - 2, Sigma = Psi)[,,1])
  sigma_sq <- diag(cov_mat) / (N - 2)
  
  # SMD
  smd <- as.vector(meandiff / sqrt(sigma_sq))  # cohen's d 
  smd <- smd * (1 - (3/((4 * (N - 2)) - 1))) # hedges g
  var_smd <- 4 / N + smd^2 / (2 * (N - 2))
  
  dat <- tibble(smd = smd, var_smd = var_smd)
  
  return(dat)
}


set.seed(202043)
# one study, 3 effect sizes 
generate_rsmd(delta = 0, k = 3, N = 70, Psi = 0.8)

# tau from Isq (Piggot 2012)
sqrt(.33* .05)
sqrt(.05)
sqrt(1.33 * .05)

# Generate meta data ------------------------------------------------------

# m is the total number of studies
# tau is the variance of the deltas
# k_mean is the mean number of effect sizes per study.
# N_mean is the mean number of participants?
# rho is the average correlation between pairs of effect sizes within a study
# nu and rho control variability of r_k across studies smaller nu more variable correlations
# covs is the design matrix
# beta is regression coefficients for X (do I need to make an intercept?)
# v_j random normal with mean 0 and variance tau sq


# generate es num id

generate_es_num <- function(dat) {
  
  dat <- dat %>%  
    group_by(study) %>%
    mutate(es_num = sequence(n())) %>%
    ungroup() %>%
    select(study, es_num, everything())
  
  return(dat)
}


# Tipton Pusto design matrix cleaned - clean_design_mat.R
load("data/design_mat.Rdata")

generate_rmeta <- function(m, tau, k_mean, N_mean, 
                           rho, nu, covs, beta,
                           return_study_params = FALSE) {
  
  # Design matrix -----------------------------------------------------------
  
  # if m > 20 repeat the design matrix 
  if (m > 20) {
    design_mat_all <- map_dfr(seq_len(ceiling(m/20)), function(x) covs)
  } else {
    design_mat_all <- covs
  }
  
  all_m <- m * 10
  
  design_mat_all <- design_mat_all %>%
    slice(1:all_m) %>%
    mutate(study = rep(1:m, each = 10)) %>% 
    generate_es_num()
  
  # Study data --------------------------------------------------------------
  
  study_data <- 
    tibble(
      k = pmin(1 + rpois(m, k_mean - 1), 10), # look at some meta analysis 
      N = pmin(20 + 2 * rpois(m, 30), 200), # distribution of sample size 
      Psi = rbeta(m, rho * nu, (1 - rho) * nu) # you have to make something up 
    ) %>%
    mutate(study = 1:m)
  
  
  # True delta --------------------------------------------------------------
  
  X <- design_mat_all %>%
    select(starts_with("X")) %>%
    as.matrix()
  
  # v_j and u_ij terms
  v_j <- rnorm(m, 0, tau)  
  study_id <- design_mat_all$study
  v_j_long <- v_j[study_id]
  #u_ij <- rnorm(m * 10, 0, omega)
  
  true_delta <- tibble(delta = as.vector(X %*% beta) + v_j_long, # is this right??
                       study = study_id) %>%
    group_by(study) %>% 
    summarize(delta = list(delta)) %>% # unnest does mini tibbles so I did this to get vector
    ungroup() 
  
  # join with the study data 
  study_data <- study_data %>%
    left_join(true_delta, by = "study") %>%
    select(-study)
  
  if (return_study_params) return(study_data)
  
  # Generate fulll meta data  -----------------------------------------------
  
  # first line runs generate_rsmd
  meta_cov_dat <- pmap_df(study_data, generate_rsmd, .id = "study") %>%
    mutate(study = as.numeric(study)) %>% # have to do study and es_num here again to join the covs
    generate_es_num() %>%
    left_join(design_mat_all, by = c("study", "es_num"))
  
  
  return(meta_cov_dat)
}


# generate meta data 
set.seed(342020)


# some of the smds are really high lol
meta_data <- 
  generate_rmeta(m = 75, 
                 tau = 0.4, 
                 k_mean = 5, 
                 N_mean = 40, 
                 rho = 0.8, 
                 nu = 39,
                 covs = design_mat,
                 beta = matrix(c(1, .1, .5, .3, .6, .7), nrow = 6))

check <- meta_data %>%
  group_by(study) %>%
  summarize(n = n())


save(meta_data, file = "data/meta_data_practice.RData")


big_meta <- 
  generate_rmeta(m = 1000, 
                 tau = 0.1, 
                 k_mean = 5, 
                 N_mean = 100, 
                 rho = 0.8, 
                 nu = 50,
                 covs = design_mat,
                 beta = matrix(c(0.3, rep(0,5)), nrow = 6),
                 return_study_params = FALSE)

mean(big_meta$var_smd)
hist(big_meta$var_smd)

# Study design features ---------------------------------------------------


# Use return_study_params to get 
# distribution of study design features

study_features <- 
  generate_rmeta(m = 1000, 
                 tau = 0.05, 
                 k_mean = 5, 
                 N_mean = 40, 
                 rho = 0.6, 
                 nu = 50,
                 covs = design_mat,
                 beta = matrix(c(0.3, rep(0,5)), nrow = 6),
                 return_study_params = TRUE)

# check means of k, N, Psi

hist(study_features$k)
hist(study_features$N)
hist(study_features$Psi)


sd(study_features$Psi) # Should be sqrt(rho (1 - rho) / nu)
sqrt(0.6 * 0.4 / 50)



# check distribution of delta

delta_data <- 
  study_features %>%
  mutate(study = 1:n()) %>%
  unnest(cols = delta)


library(lme4)
lmer(delta ~ (1 | study), data = delta_data)
# Estimated SDs should be very close to tau and omega


# True parameter recovery -------------------------------------------------



library(clubSandwich)
library(metafor)

# Set nu = 4000 so that all Psi = rho

meta_data <- 
  generate_rmeta(m = 1000, 
                 tau = 0.05, 
                 k_mean = 5, 
                 N_mean = 40, 
                 rho = 0.6, 
                 nu = 4000,
                 covs = design_mat,
                 beta = matrix(c(1, .1, .5, .3, .6, .7), nrow = 6))

# save(meta_data, file = "data/meta_data_practice.Rdata")

# Use clubSandwich::impute_covariance_matrix with true rho.
V_mat <- impute_covariance_matrix(vi = meta_data$var_smd, 
                                  cluster = meta_data$study, 
                                  r = 0.6)

