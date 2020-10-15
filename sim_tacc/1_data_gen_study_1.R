#------------------------------------------------------
# Data Generating Model
#------------------------------------------------------



generate_rsmd <- function(delta, k, N, Sigma) {
  
  # make sure delta is a vector
  delta_vec <- rep(delta, length.out = k)
  
  # create Sigma matrix assuming equicorrelation
  if (!is.matrix(Sigma)) Sigma <- Sigma + diag(1 - Sigma, nrow = k) # cor matrix for 1 study
  
  # generate numerator of SMD 
  mean_diff <- rmvnorm(n = 1, mean = delta_vec, sigma = (4/N) * Sigma) 
  
  # covariance 
  cov_mat <- as.matrix(rWishart(n = 1, df = N - 2, Sigma = Sigma)[,,1])
  sigma_sq <- diag(cov_mat) / (N - 2)
  
  # SMD
  d <- as.vector(mean_diff / sqrt(sigma_sq))  # cohen's d 
  J <- (1 - (3/((4 * (N - 2)) - 1)))
  g <- d * J # hedges g
  var_g <- J^2 * (4 / N + d^2 / (2 * (N - 2)))
  
  dat <- tibble(g = g, var_g = var_g)
  
  return(dat)
}


generate_es_num <- function(dat) {
  
  dat <- dat %>%  
    group_by(study) %>%
    mutate(es_num = sequence(n())) %>%
    ungroup() %>%
    select(study, es_num, everything())
  
  return(dat)
}




generate_rmeta <- function(m, 
                           tau, 
                           rho, 
                           covs, 
                           beta_type = "A",
                           k_mean = 4,
                           N_mean = 30,
                           nu = 50,
                           return_study_params = FALSE) {
  
  
  
  # beta --------------------------------------------------------------------
  
  if(beta_type == "A") {
    beta <- c(.3, 0, 0, 0, 0, 0)
  } else if(beta_type == "B1"){
    beta <- c(.3, .1, 0, 0, 0, 0)
  } else if(beta_type == "B5"){
    beta <- c(.3, .5, 0, 0, 0, 0)
  } else if(beta_type == "C1"){
    beta <- c(.3, 0, .1, 0, 0, 0)
  } else if(beta_type == "C5"){
    beta <- c(.3, 0, .5, 0, 0, 0)
  } else if(beta_type == "D1"){
    beta <- c(.3, 0, 0, .1, 0, 0)
  } else if(beta_type == "D5"){
    beta <- c(.3, 0, 0, .5, 0, 0)
  } else if(beta_type == "E1"){
    beta <- c(.3, 0, 0, 0, .1, 0)
  } else if(beta_type == "E5"){
    beta <- c(.3, 0, 0, 0, .5, 0)
  } else if(beta_type == "F1"){
    beta <- c(.3, 0, 0, 0, 0, .1)
  } else if(beta_type == "F5"){
    beta <- c(.3, 0, 0, 0, 0, .5)
  } 
  
  beta <- matrix(beta, nrow = 6)
  
  
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
      k = pmin(1 + rpois(m, k_mean), 10), # look at some meta analysis 
      N = pmin(20 + 2 * rpois(m, N_mean), 200), # distribution of sample size 
      Sigma = rbeta(m, rho * nu, (1 - rho) * nu) # you have to make something up 
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
  
  # Generate full meta data  -----------------------------------------------
  
  # first line runs generate_rsmd
  meta_cov_dat <- 
    pmap_df(study_data, generate_rsmd, .id = "study") %>%
    mutate(study = as.numeric(study)) %>% # have to do study and es_num here again to join the covs
    generate_es_num() %>%
    left_join(design_mat_all, by = c("study", "es_num"))
  
  
  return(meta_cov_dat)
}

