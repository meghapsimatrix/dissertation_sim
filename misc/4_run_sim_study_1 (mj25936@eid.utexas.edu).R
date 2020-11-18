library(dplyr)
library(purrr)
library(mvtnorm)
library(robumeta)
library(clubSandwich)
library(tidyr)

# Tipton Pusto design matrix cleaned - clean_design_mat.R
load("data/design_mat.Rdata")
load("data/to_test.Rdata")

#-----------------------------------------------------------
# Source the functions
#-----------------------------------------------------------

source("1_data_gen_study_1.R")
source("2_estimation_study_1.R")
source("3_performance_criteria.R")

#-----------------------------------------------------------
# Simulation Driver - should return a data.frame or tibble
#-----------------------------------------------------------

# JEP: It might be easier to pare down the test_dat outside of this function,
#      as part of the design for the power sims. Then you could nest() the 
#      test_dat and pass it as an argument.

# instead of passing to_test - simulation design as a parameter dataset 
# make a little data frame with just beta types 
# Cross it with full to_test filter to just the tests that you want 
# nest it one row per beta type 

# put R in params instead of to_test

run_sim <- function(iterations, m, tau, rho, beta_type, R, 
                    full_form = "g ~ X1 + X2 + X3 + X4 + X5", 
                    design_matrix = design_mat, 
                    test_dat = to_test, 
                    seed = NULL) {
  
    
  if (!is.null(seed)) set.seed(seed)
  
  results <-
    map_dfr(1:iterations, {
        
      # generate data ------------------------------------------------------------
      meta_data <- generate_rmeta(m = m, tau = tau, rho = rho, beta_type = beta_type)
      
      # Fit full model on data --------------------------------------------------
      full_formula <- full_form
      full_model <- robu(as.formula(full_formula), 
                         studynum = study, 
                         var.eff.size = var_g,
                         small = FALSE,
                         data = meta_data)
      
      # get cov matrices --------------------------------------------------------
      cov_mat_cr1 <- vcovCR(full_model, type = "CR1")
      cov_mat_cr2 <- vcovCR(full_model, type = "CR2")
      
      # get naive and htz -------------------------------------------------------
      names(test_dat$indices_test) <- test_dat$cov_test
      naive_res <- Wald_test(full_model, 
                             constraints = constrain_zero(test_dat$indices_test),
                             vcov = cov_mat_cr1,
                             test = "Naive-F", 
                             tidy = TRUE) %>%
        select(`Naive-F` = p_val)
      
      
      htz_res <- Wald_test(full_model, 
                           constraints = constrain_zero(test_dat$indices_test),
                           vcov = cov_mat_cr2,
                           test = "HTZ", 
                           tidy = TRUE) %>%
        select(HTZ = p_val)
      

      # cwb ---------------------------------------------------------------------
      # put R in params
      cwb_params <- test_dat %>%
        select(null_model, indices_test, R)
        select(null_model, indices_test, beta_seed)

      
      boot_res <- pmap_dfr(cwb_params, cwb)
      
      res <- 
        bind_cols(naive_res, htz_res, boot_res) %>%
        bind_cols(test_dat %>% select(cov_test, contrasts)) %>%
        gather(test, p_val, -c(cov_test, contrasts))
      
    }) 
  
  calc_performance(res)
}

# demonstrate the simulation driver

#-------------------------------------
# Experimental Design
#-------------------------------------
source_obj <- ls()

set.seed(20150316) # change this seed value!

# now express the simulation parameters as vectors/lists

design_factors <- list(
  m = c(10, 20, 40, 80), 
  tau = c(0.1, 0.3),
  rho = c(0.5, 0.8),
  beta_type = c("A", "B1", "B5", "C1", "C5", "D1", "D5", "E1", "E5", "F1", "F5")
)

# combine into a design set

params <-
  cross_df(design_factors) %>%
  mutate(
    iterations = 4000, # change this to how many ever iterations
    seed = round(runif(1) * 2^30) + 1:n()
  )


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

save(params, results, session_info, run_date, file = "simulation_results_study_1.Rdata")