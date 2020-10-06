#------------------------------------------------------
# Set development values for simulation parameters
#------------------------------------------------------

# What are your model parameters?
# What are your design parameters?

library(dplyr)
library(purrr)
library(mvtnorm)


#------------------------------------------------------
# Data Generating Model
#------------------------------------------------------

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
  d <- as.vector(meandiff / sqrt(sigma_sq))  # cohen's d 
  J <- (1 - (3/((4 * (N - 2)) - 1)))
  g <- d * (1 - (3/((4 * (N - 2)) - 1))) # hedges g
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


# Tipton Pusto design matrix cleaned - clean_design_mat.R
load("data/design_mat.Rdata")

generate_rmeta <- function(m, tau, 
                           rho, covs, beta_type,
                           return_study_params = FALSE) {
  
  
  # mean es num and N -------------------------------------------------------
  
  k_mean <- 4
  N_mean <- 30
  nu <- 50
  
  
  # beta --------------------------------------------------------------------
  
  if(beta_type == "A") {
    beta <- c(.3, 0, 0, 0, 0, 0, 0, 0)
  } else if(beta_type == "B1"){
    beta <- c(.3, .1, 0, 0, 0, 0, 0, 0)
  } else if(beta_type == "B5"){
    beta <- c(.3, .5, 0, 0, 0, 0, 0, 0)
  } else if(beta_type == "C1"){
    beta <- c(.3, 0, .1, 0, 0, 0, 0, 0)
  } else if(beta_type == "C5"){
    beta <- c(.3, 0, .5, 0, 0, 0, 0, 0)
  } else if(beta_type == "D1"){
    beta <- c(.3, 0, 0, .1, 0, 0, 0, 0)
  } else if(beta_type == "D5"){
    beta <- c(.3, 0, 0, .5, 0, 0, 0, 0)
  } else if(beta_type == "E1"){
    beta <- c(.3, 0, 0, 0, .1, 0, 0, 0)
  } else if(beta_type == "E5"){
    beta <- c(.3, 0, 0, 0, .5, 0, 0, 0)
  } else if(beta_type == "F1"){
    beta <- c(.3, 0, 0, 0, 0, .1, 0, 0)
  } else if(beta_type == "F5"){
    beta <- c(.3, 0, 0, 0, 0, .5, 0, 0)
  } 
  
  beta <- matrix(beta, nrow = 8)
  
  
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



# Test the data-generating model - How can you verify that it is correct?


#------------------------------------------------------
# Model-fitting/estimation/testing functions
#------------------------------------------------------


estimate_wald <- function(model, indices_test, cov_mat, test){
  
  res <- Wald_test(model, 
            constraints = constrain_zero(indices_test), 
            vcov = cov_mat,
            test = test) %>%
    select(test, p_val)
  
  return(res)
}




# Test the estimation function

#------------------------------------------------------
# Calculate performance measures
# (For some simulations, it may make more sense
# to do this as part of the simulation driver.)
#------------------------------------------------------

calc_performance <- function(results, alpha) {
  
  performance_measures  <- results %>%
    filter(!is.na(p_val)) %>%
    group_by(test) %>%
    summarize(K = n(),
              rej_rate = mean(p_val < alpha),
              mcse = sqrt((rej_rate * (1 - rej_rate))/K))

  return(performance_measures)
}

# Check performance calculations

#-----------------------------------------------------------
# Simulation Driver - should return a data.frame or tibble
#-----------------------------------------------------------

run_sim <- function(iterations, model_params, design_params, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)

  results <-
    rerun(iterations, {
      
      dat <- generate_rmeta(model_params)
      naive_res <- map_dfr(to_test$indices_test, estimate_wald, model = full_model, cov_mat = cov_mat_cr1, test = "Naive-F")
      htz_res <- map_dfr(to_test$indices_test, estimate_wald, model = full_model, cov_mat = cov_mat_cr2, test = "HTZ")
      results <- bind_rows(naive_res, htz_res, boot_res)
      
    }) %>%
    bind_rows()

  calc_performance(results, alpha)
}

# demonstrate the simulation driver

#-------------------------------------
# Experimental Design
#-------------------------------------
set.seed(20150316) # change this seed value!

# now express the simulation parameters as vectors/lists

beta <- list(c(.3, 0, 0, 0, 0, 0, 0, 0),
             c(.3, .1, 0, 0, 0, 0, 0, 0),
             c(.3, .5, 0, 0, 0, 0, 0, 0),
             c(.3, 0, .1, 0, 0, 0, 0, 0),
             c(.3, 0, .5, 0, 0, 0, 0, 0), 
             c(.3, 0, 0, .1, 0, 0, 0, 0), 
             c(.3, 0, 0, .5, 0, 0, 0, 0),
             c(.3, 0, 0, 0, .1, 0, 0, 0), 
             c(.3, 0, 0, 0, .5, 0, 0, 0),
             c(.3, 0, 0, 0, 0, .1, 0, 0), 
             c(.3, 0, 0, 0, 0, .5, 0, 0), 
             c(.3, 0, 0, 0, 0, 0, .1, 0),
             c(.3, 0, 0, 0, 0, 0, .5, 0), 
             c(.3, 0, 0, 0, 0, 0, 0, .1),
             c(.3, 0, 0, 0, 0, 0, 0, .5))

design_factors <- list(
  m = c(10, 20, 40, 80), 
  tau = c(0.1, 0.3),
  rho = c(0.5, 0.8),
  beta_type = c("A", "B1", "B5", "C1", "C5", "D1", "D5", "E1", "E5", "F1", "F5", "G1", "G5", "H1", "H5")
  )

# combine into a design set

params <-
  cross_df(design_factors) %>%
  mutate(
    iterations = 1000, # change this to how many ever iterations
    seed = round(runif(1) * 2^30) + 1:n()
  )

# All look right?
lengths(design_factors)
nrow(params)
head(params)


#--------------------------------------------------------
# run simulations in parallel - mdply workflow
#--------------------------------------------------------

library(Pusto)

cluster <- start_parallel(source_obj = source_obj, register = TRUE)

system.time(results <- plyr::mdply(params, .fun = run_sim, .parallel = TRUE))

stopCluster(cluster)


#--------------------------------------------------------
# Save results and details
#--------------------------------------------------------

session_info <- sessionInfo()
run_date <- date()

save(params, results, session_info, run_date, file = "simulation_results.Rdata")