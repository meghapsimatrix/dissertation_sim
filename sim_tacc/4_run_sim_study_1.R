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

run_sim <- function(iterations, m, tau, rho, beta_type, design_matrix, test_dat = to_test, R, boot_seed, seed = NULL) {
  

  # non zero betas only for power -------------------------------------------
      
  if (beta_type %in% c("B1", "B5")){
        
      test_dat <- test_dat %>%
        filter(str_detect(cov_test, "X1"))
        
  } else if(beta_type %in% c("C1", "C5")){
        
      test_dat <- test_dat %>%
        filter(str_detect(cov_test, "X2"))
      
  } else if(beta_type %in% c("D1", "D5")){
        
    test_dat <- test_dat %>%
      filter(str_detect(cov_test, "X3"))
        
  } else if(beta_type %in% c("E1", "E5")){
        
    test_dat <- test_dat %>%
      filter(str_detect(cov_test, "X4"))
        
  } else if(beta_type %in% c("F1", "F5")){
        
    test_dat <- test_dat %>%
      filter(str_detect(cov_test, "X5"))
        
  } else if (beta_type == "A" ){
        
    test_dat <- test_dat 
        
  } 
      
    
  if (!is.null(seed)) set.seed(seed)
  

  
  results <-
    map_dfr(1:iterations, {
        
      # generate data ------------------------------------------------------------
      meta_data <- generate_rmeta(m = m, tau = tau, rho = rho, beta_type = beta_type)
      
      # Fit full model on data --------------------------------------------------
      full_formula <- "g ~ X1 + X2 + X3 + X4 + X5"
      full_model <- robu(as.formula(full_formula), 
                         studynum = study, 
                         var.eff.size = var_g,
                         small = FALSE,
                         data = meta_data)
      
      # get cov matrices --------------------------------------------------------
      cov_mat_cr1 <- vcovCR(full_model, type = "CR1")
      cov_mat_cr2 <- vcovCR(full_model, type = "CR2")
      
      # get naive and htz -------------------------------------------------------
      
      # JEP: This can take advantage of the mapping features built-in to Wald_test()
      #      See comments in tests/2_est_perf_test.R
      
      naive_res <- map_dfr(test_dat$indices_test,
                           estimate_wald, 
                           model = full_model, 
                           cov_mat = cov_mat_cr1, 
                           test = "Naive-F")
      
      
      htz_res <-  map_dfr(test_dat$indices_test,
                          estimate_wald, 
                          model = full_model, 
                          cov_mat = cov_mat_cr2, 
                          test = "HTZ")

      # cwb ---------------------------------------------------------------------
      # JEP: See suggestions in 2_estimation_study_1.R regarding fitting the null
      #      Model inside of cwb().
      
      null_mods <- map(test_dat$null_model, fit_mod)
      
      cwb_params <- test_dat %>%
        mutate(null_mod = null_mods) %>%
        select(null_mod, indices_test)
      
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