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

   dat$es_num <- sequence(rle(dat$study)$lengths)

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
    beta <- c(0, 0, 0, 0, 0, 0)
  } else if(beta_type == "B5"){
    beta <- c(0, .5, 0, 0, 0, 0)
  }

  
  beta <- matrix(beta, nrow = 6)
  
  
  # Study data --------------------------------------------------------------
  
  study_data <- 
    tibble(
      k = pmin(1 + rpois(m, k_mean), 10), # look at some meta analysis 
      N = pmin(20 + 2 * rpois(m, N_mean), 200), # distribution of sample size 
      Sigma = rbeta(m, rho * nu, (1 - rho) * nu) # you have to make something up 
    ) %>%
    mutate(study = 1:m)
  
  
  # Design matrix -----------------------------------------------------------

  cat <- LETTERS[1:5]
  X1 <- sample(cat, m, replace = TRUE, prob = c(0.2, 0.2, 0.2, 0.2, 0.2))
  X2 <- sample(cat, sum(study_data$k), replace = TRUE, prob = c(0.2, 0.2, 0.2, 0.2, 0.2))
  
  design_mat_all <- tibble(X = 1, 
                           X1 = rep(X1, study_data$k), 
                           X2 = X2) %>%
    mutate(study = rep(1:m, study_data$k)) %>% 
    generate_es_num() %>%
    dplyr::select(-X2)
  
  design_mat_all <- dummy_cols(design_mat_all, select_columns = "X1",
                               remove_selected_columns = TRUE)
  
  # True delta --------------------------------------------------------------
  
  X <- design_mat_all %>%
    dplyr::select(starts_with("X")) %>%
    as.matrix()
  
  
  # v_j 
  v_j <- rnorm(m, 0, tau)  
  study_id <- design_mat_all$study
  v_j_long <- v_j[study_id]
  ##u_ij <- rnorm(m * 10, 0, omega)
  
  true_delta <- tibble(delta = as.vector(X %*% beta) + v_j_long, # is this right??
                       study = study_id) %>%
    group_by(study) %>% 
    summarize(delta = list(delta), .groups = "drop") %>% # unnest does mini tibbles so I did this to get vector
    ungroup() 
  
  # join with the study data 
  study_data <- study_data %>%
    dplyr::left_join(true_delta, by = "study") %>%
    dplyr::select(-study)
  
  if (return_study_params) return(study_data)
  
  # Generate full meta data  -----------------------------------------------
  
  # first line runs generate_rsmd
  meta_cov_dat <- 
    pmap_df(study_data, generate_rsmd, .id = "study") %>%
    mutate(study = as.numeric(study)) %>% # have to do study and es_num here again to join the covs
    generate_es_num() %>%
    dplyr::left_join(design_mat_all, by = c("study", "es_num"))
  
  return(meta_cov_dat)
}
