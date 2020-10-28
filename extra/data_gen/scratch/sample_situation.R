sample(1:10, 75, replace = TRUE)

m <- 75
rho <- .6
nu <- 50
covs <- design_mat
tau <- 0.05
beta <- matrix(c(1, .1, .5, .3, .6, .7), nrow = 6)


study_data <- 
  tibble(
    k = sample(1:10, m, replace = TRUE), # 1:10
    N = sample(32:130, m, replace = TRUE), # distribution of sample size 
    Psi = rbeta(m, rho * nu, (1 - rho) * nu) # you have to make something up 
  ) %>%
  mutate(study = 1:m)


design_mat_all <- map_dfr(seq_len(ceiling(m/20)), function(x) covs)


design_mat_all <- design_mat_all %>%
  slice(1:all_m) %>%
  mutate(study = rep(1:m, each = 10)) %>% 
  generate_es_num()

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


meta_cov_dat <- pmap_df(study_data, generate_rsmd, .id = "study") 




delta <- unlist(study_data$delta[1])
k <- study_data$k[1]
N <- study_data$N[1]
Psi <- study_data$Psi[1]


# make sure delta is a vector
delta_vec <- rep(delta, length.out = k)

# create Psi matrix assuming equicorrelation
if (!is.matrix(Psi)) Psi <- Psi + diag(1 - Psi, nrow = k) # cor matrix for 1 study

# generate numerator of SMD 
meandiff <- rmvnorm(n = 1, mean = delta_vec, sigma = (4/N) * Psi) 

# covariance 
cov_mat <- rWishart(n = 1, df = N - 2, Sigma = Psi)[,,1]
sigma_sq <- diag(cov_mat) / (N - 2)

# SMD
smd <- as.vector(meandiff / sqrt(sigma_sq))
v <- sqrt(4 / N + smd^2 / (2 * (N - 2)))

dat <- tibble(smd = smd, v = v)
