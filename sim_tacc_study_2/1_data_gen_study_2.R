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
                           cov_type,
                           cat_num,
                           beta_type = "A",
                           k_mean = 4,
                           N_mean = 30,
                           nu = 50,
                           return_study_params = FALSE) {
  
  
  
  # beta --------------------------------------------------------------------
  
  # JAMES PLEASE CHECK

  if (beta_type == "A") {
    beta <- c(.3, 0, rep(0, cat_num - 2))
  } else if(beta_type == "B5"){
    beta <- c(.3, .5, rep(0, cat_num - 2))
  }

  
  beta <- matrix(beta, nrow = cat_num)
  
  
  # Study data --------------------------------------------------------------
  
  study_data <- 
    tibble(
      k = pmin(1 + rpois(m, k_mean), 10), # look at some meta analysis 
      N = pmin(20 + 2 * rpois(m, N_mean), 200), # distribution of sample size 
      Sigma = rbeta(m, rho * nu, (1 - rho) * nu) # you have to make something up 
    ) %>%
    mutate(study = 1:m)
  
  
  # Design matrix -----------------------------------------------------------

  # JAMES PLEASE CHECK
  
  cat <- LETTERS[1:cat_num]
  min_times <- 2
  
  if(cov_type == "between") {
  
    cat_var <- c(rep(cat, each = min_times), sample(cat, size = m - min_times * cat_num, replace = TRUE))
    X1 <- sample(cat_var)
    
    design_mat_all <- tibble(X = 1, 
                             X1 = rep(X1, study_data$k)) 
    
  } else if(cov_type == "within"){
    
    cat_var <- c(rep(cat, each = min_times), sample(cat, size = sum(study_data$k) - min_times * cat_num, replace = TRUE))
    X1 <- sample(cat_var)
    
    design_mat_all <- tibble(X = 1, 
                             X1 = X1) 
    
  }
  
  design_mat_all <- design_mat_all %>%
    mutate(study = rep(1:m, study_data$k)) %>% 
    generate_es_num()
  
  design_mat_all <- dummy_cols(design_mat_all, select_columns = "X1",
                               remove_selected_columns = TRUE) %>%
    select(-X1_A)
  
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
