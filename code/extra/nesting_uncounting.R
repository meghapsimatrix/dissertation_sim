library(tidyverse)

m = 50
mu = 0.2 
tau = 0.05 
k_mean = 5 
N_mean = 40 
rho = 0.6 
nu = 39

set.seed(202045)

study_data <- 
  tibble(
    delta = rnorm(m, mean = mu, sd = tau),
    k = pmin(2 + rpois(m, k_mean - 2), 10), # is max 10, I don't understand this part at all :D 
    N = 20 + 2 * rpois(m, (N_mean - 20) / 2),
    Psi = rbeta(m, rho * nu, (1 - rho) * nu)
  )

study_data_all <- study_data %>%
  uncount(k) %>%
  mutate(study = rep(1:m, times = study_data$k))


# Tipton Pusto design matrix cleaned - clean_design_mat.R
load("data/design_mat.Rdata")

covs <- design_mat

design_mat_all <- map_dfr(seq_len(ceiling(m/20)), function(x) covs)

# m_max <- 20 * ceiling(m / 2) 

all_m <- m * 10

design_mat_all <- design_mat_all %>%
  slice(1:all_m) %>%
  mutate(study = rep(1:m, each = 10)) %>%
  group_by(study) %>%
  mutate(es_num = sequence(n())) %>%
  ungroup()

meta_all <- left_join(study_data_all, design_mat_all, by = c("study", "es_num"))


