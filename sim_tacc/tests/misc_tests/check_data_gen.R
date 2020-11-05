library(dplyr)
library(purrr)
library(mvtnorm)
library(robumeta)
library(clubSandwich)
library(tidyr)


load("data/design_mat.Rdata")

#-----------------------------------------------------------
# Source the functions
#-----------------------------------------------------------

source("1_data_gen_study_1.R")


covs <- design_mat
m <- 80

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

k_mean <- 4
N_mean <- 30
rho <- .6
nu <- 50
tau <- .1

study_data <- 
  tibble(
    k = pmin(1 + rpois(m, k_mean), 10), # look at some meta analysis 
    N = pmin(20 + 2 * rpois(m, N_mean), 200), # distribution of sample size 
    Sigma = rbeta(m, rho * nu, (1 - rho) * nu) # you have to make something up 
  ) %>%
  mutate(study = 1:m)

beta <- c(.3, 0, 0, 0, 0, .5)
beta <- matrix(beta, nrow = 6)

X <- design_mat_all %>%
  dplyr::select(starts_with("X")) %>%
  as.matrix()

# v_j and u_ij terms
v_j <- rnorm(m, 0, tau)  
study_id <- design_mat_all$study
v_j_long <- v_j[study_id]
#u_ij <- rnorm(m * 10, 0, omega)


true_delta <- tibble(delta = as.vector(X %*% beta) + v_j_long, # is this right??
                     study = study_id)

delta <- split(true_delta$delta, true_delta$study)
k <- split(study_data$k, study_data$study)
N <- split(study_data$N, study_data$study)
Sigma <- split(study_data$Sigma, study_data$study)



meta_cov_dat <- 
  pmap_dfr(list(delta, k, N, Sigma), generate_rsmd, .id = "study") %>%
  mutate(study = as.numeric(study)) %>% # have to do study and es_num here again to join the covs
  generate_es_num() %>%
  dplyr::left_join(design_mat_all, by = c("study", "es_num"))

glimpse(meta_cov_dat)
