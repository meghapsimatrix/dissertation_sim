# Generate meta data ------------------------------------------------------

# m is the total number of studies
# mu is the mean of the deltas
# tau is the variance of the deltas
# k_mean is the mean number of effect sizes per study.
# N_mean is the mean number of participants?
# rho is theaverage correlation between pairs of effect sizes within a study
# nu and rho control variability of r_k across studies smaller nu more variable correlations
# covs is the design matrix


# Tipton Pusto design matrix cleaned - clean_design_mat.R
load("data/design_mat.Rdata")

generate_rmeta <- function(m, mu, tau, k_mean, N_mean, rho, nu, covs) {
  
  # design matrix
  # if m > 20 repeat the design matrix 
  if(m > 20) {
    design_mat_all <- map_dfr(seq_len(ceiling(m/20)), function(x) covs)
  } else{
    design_mat_all <- covs
  }
  
  all_m <- m * 10
  
  design_mat_all <- design_mat_all %>%
    slice(1:all_m) %>%
    mutate(study = rep(1:m, each = 10)) %>%
    group_by(study) %>%
    mutate(es_num = sequence(n())) %>%
    ungroup()
  
  
  # study data 
  study_data <- 
    tibble(
      delta = rnorm(m, mean = mu, sd = tau),
      k = pmin(2 + rpois(m, k_mean - 2), 10), 
      N = 20 + 2 * rpois(m, (N_mean - 20) / 2),
      Psi = rbeta(m, rho * nu, (1 - rho) * nu)
    )
  
  study_data_all <- study_data %>%
    uncount(k) %>%
    mutate(study = rep(1:m, times = study_data$k))
  
  
  # meta data 
  meta_dat <- pmap_df(study_data, generate_rsmd, .id = "study") %>%
    mutate(study = as.numeric(study)) %>%
    group_by(study) %>%
    mutate(es_num = sequence(n())) %>%
    ungroup()
  

  meta_cov_dat <- left_join(meta_dat, design_mat_all, by = c("study", "es_num")) %>%
    select(study, es_num, everything())
  
  return(meta_cov_dat)
}


# generate meta data 
set.seed(342020)

meta_data <- 
  generate_rmeta(m = 75, 
                 mu = 0.2, 
                 tau = 0.05, 
                 k_mean = 5, 
                 N_mean = 40, 
                 rho = 0.6, 
                 nu = 39,
                 covs = design_mat)

check <- meta_data %>%
  group_by(study) %>%
  summarize(n = n())
