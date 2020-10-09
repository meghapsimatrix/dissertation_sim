library(dplyr)
library(purrr)
library(mvtnorm)
library(robumeta)
library(clubSandwich)

# Tipton Pusto design matrix cleaned - clean_design_mat.R
load("data/design_mat.Rdata")

#-----------------------------------------------------------
# Source the functions
#-----------------------------------------------------------

source("1_data_gen_study_1.R")
source("3_performance_criteria.R")

#-----------------------------------------------------------
# Simulation Driver - should return a data.frame or tibble
#-----------------------------------------------------------

run_sim <- function(iterations, model_params, design_params, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  results <-
    rerun(iterations, {
      

      # generate data ------------------------------------------------------------
      meta_data <- generate_rmeta(model_params)
      

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
      
      naive_res <- map_dfr(to_test$indices_test,
                           estimate_wald, 
                           model = full_model, 
                           cov_mat = cov_mat_cr1, 
                           test = "Naive-F")
      
      
      htz_res <-  map_dfr(to_test$indices_test,
                          estimate_wald, 
                          model = full_model, 
                          cov_mat = cov_mat_cr2, 
                          test = "HTZ")
      

      # cwb ---------------------------------------------------------------------

      null_mods <- map(to_test$null_model, fit_mod)
      
      cwb_params <- to_test %>%
        mutate(null_mod = null_mods) %>%
        select(null_mod, indices_test)
      
      boot_res <- pmap_dfr(cwb_params, cwb)
      
      cwb_res <- boot_res %>%
        filter(test == "CWB") 
      
      cwb_a_res <- boot_res %>%
        filter(test == "CWB Adjusted") 
      
      res <- bind_cols(naive_res, htz_res, cwb_res, cwb_a_res) %>%
        bind_cols(to_test %>% select(cov_test)) %>%
        gather(test, p_val, -cov_test)
      
    }) %>%
    bind_rows()
  
  calc_performance(results, alpha)
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