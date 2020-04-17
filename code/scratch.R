library(tidyverse)
library(mvtnorm)

load("data/design_mat.Rdata")


m <- 75
tau <- 0.05
k_mean <- 5
N_mean <- 40
rho <- 0.6
nu <- 39
covs <- design_mat
beta <- matrix(c(1, .1, .5, .3, .6, .7), nrow = 6)
covs <- design_mat

generate_es_num <- function(dat){
  
  dat <- dat %>%  
    group_by(study) %>%
    mutate(es_num = sequence(n())) %>%
    ungroup() %>%
    select(study, es_num, everything())
  
  return(dat)
}  

load("data/design_mat.Rdata")
design_mat_all <- map_dfr(seq_len(ceiling(m/20)), function(x) covs)

all_m <- m * 10

design_mat_all <- design_mat_all %>%
  slice(1:all_m) %>%
  mutate(study = rep(1:m, each = 10)) %>% 
  generate_es_num()

# Study data --------------------------------------------------------------

study_data <- 
  tibble(
    k = pmin(2 + rpois(m, k_mean - 2), 10), 
    N = 20 + 2 * rpois(m, (N_mean - 20) / 2),
    Psi = rbeta(m, rho * nu, (1 - rho) * nu)
  ) %>%
  mutate(study = 1:m)



X <- design_mat_all %>%
  select(starts_with("X")) %>%
  as.matrix()

v_j <- rnorm(m, 0, tau)  
study_id <- rep(1:m, each = 10)
vj_long <- v_j[study_id]

true_delta <- tibble(delta = as.vector(X %*% beta) + v_j, # is this right??
                     study = study_id) %>%
  group_by(study) %>% 
  summarize(delta = list(delta)) %>% # unnest does mini tibbles so I did this to get vector
  ungroup() 


study_data <- study_data %>%
  left_join(true_delta, by = "study") %>%
  select(-study)

# first line runs generate_rsmd
meta_cov_dat <- pmap_df(study_data, generate_rsmd, .id = "study") %>%
  mutate(study = as.numeric(study)) %>% # have to do study and es_num here again to join the covs
  generate_es_num() %>%
  left_join(design_mat_all)


