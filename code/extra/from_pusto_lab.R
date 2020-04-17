# Study data --------------------------------------------------------------

study_data <- 
  tibble(
    k = pmin(2 + rpois(m, k_mean - 2), 10), 
    N = 20 + 2 * rpois(m, (N_mean - 20) / 2),
    Psi = rbeta(m, rho * nu, (1 - rho) * nu)
  ) %>%
  mutate(study = 1:m)

# uncount it 
study_data_all <- study_data %>%   # have to do study and es num here to match with design mat
  uncount(k) %>% # uncount 
  generate_es_num()

# join the design matrix to match the study design
meta_all <- left_join(study_data_all, design_mat_all, by = c("study", "es_num")) %>%
  select(study, es_num, everything())


X <- meta_all %>%
  select(starts_with("X")) %>%
  as.matrix()


v_j <- rnorm(nrow(meta_all), 0, tau)  

true_delta <- tibble(delta = as.vector(X %*% beta) + v_j, # is this right??
                     study = meta_all$study) %>%
  group_by(study) %>% 
  summarize(delta = list(delta)) %>% # unnest does mini tibbles so I did this to get vector
  ungroup() 

study_data <- 
    data.frame(
      J = 2 + rpois(m, k_mean - 2),
      N = 20 + 2 * rpois(m, (N_mean - 20) / 2),
      Psi = rbeta(m, rho * nu, (1 - rho) * nu),
      delta_k = rnorm(m, mean = X %*% beta, sd = tau)  # what?
    )
  
  X <- as.data.frame(X) %>%
    mutate(study = as.character(1:K))
  
  pmap_df(study_data, r_SMDs_within, omega = omega, .id = "study") %>%
    left_join(X)
}