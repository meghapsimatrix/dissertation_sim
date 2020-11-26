library(dplyr)
library(purrr)
library(mvtnorm)
library(clubSandwich)
library(tidyr)
library(stringr)
library(tibble)

# batch <- commandArgs()  # command line

# Tipton Pusto design matrix cleaned - clean_design_mat.R
load("data/design_mat.Rdata")
load("data/to_test.RData")

#-----------------------------------------------------------
# Source the functions
#-----------------------------------------------------------

source("1_data_gen_study_1.R")
source("2_estimation_study_1.R")
source("3_performance_criteria.R")

#-----------------------------------------------------------
# Simulation Driver - should return a data.frame or tibble
#-----------------------------------------------------------

run_sim <- function(iterations, 
                    m, 
                    tau, 
                    rho, 
                    beta_type, 
                    batch,
                    R,
                    indices_test = 2:5,
                    full_form = "X1_B + X1_C + X1_D + X1_E", 
                    seed = NULL) {
  
  require(dplyr)
  require(purrr)
  require(mvtnorm)
  require(clubSandwich)
  require(tidyr)
  require(stringr)
  require(tibble)
  
  if (!is.null(seed)) set.seed(seed)

  
  results <-
    rerun(iterations, {
        
      # generate data ------------------------------------------------------------
      meta_data <- generate_rmeta(m = m, 
                                  tau = tau, 
                                  rho = rho,
                                  covs = design_matrix,
                                  beta_type = beta_type)
      
      # Fit full model on data --------------------------------------------------
      y <- meta_data$g
      v <- meta_data$var_g
      X <- model.matrix(as.formula(paste("g ~ ", full_form)), data = meta_data)
      cluster <- meta_data$study
      
      full_model <- robu_handmade(X = X, y = y, v = v, cluster = cluster)
      
      # get cov matrices --------------------------------------------------------
      cov_mat_cr1 <- vcovCR(full_model, type = "CR1", cluster = cluster)
      cov_mat_cr2 <- vcovCR(full_model, type = "CR2", cluster = cluster)
      
      # get naive and htz -------------------------------------------------------
      naive_res <- Wald_test(full_model, 
                             constraints = constrain_zero(indices_test),
                             vcov = cov_mat_cr1,
                             test = "Naive-F", 
                             tidy = TRUE) %>%
        dplyr::select(`Naive-F` = p_val)
      
      
      htz_res <- Wald_test(full_model, 
                           constraints = constrain_zero(indices_test),
                           vcov = cov_mat_cr2,
                           test = "HTZ", 
                           tidy = TRUE) %>%
        mutate(p_val = ifelse(Fstat < 0, 1, p_val)) %>%  # added this to fix the htz p val issue with negative F stat
        dplyr::select(HTZ = p_val)
      

      # cwb ---------------------------------------------------------------------
      system.time(boot_res <- cwb(null_model = "g ~ 1", 
                             indices_test = indices_test, 
                             full_form = full_form, 
                             R = R, dat = meta_data))
      
      res <- 
        bind_cols(naive_res, htz_res, boot_res) %>%
        bind_cols(test_dat %>% dplyr::select(cov_test)) %>%
        gather(test, p_val, -c(cov_test))
      
    })  %>%
    bind_rows()
  
  calc_performance(results)
}

# demonstrate the simulation driver

#-------------------------------------
# Experimental Design
#-------------------------------------

# include design matrix, exclude to_test

set.seed(20201108) # change this seed value!

# now express the simulation parameters as vectors/lists

design_factors <- list(
  m = c(10, 20, 40, 80), 
  tau = c(0.1, 0.3),
  rho = c(0.5, 0.8),
  R = 399,
  beta_type = c("A", "B5"),
  batch = 1:48
)

# combine into a design set
params <-
  cross_df(design_factors) %>%
  mutate(
    iterations = 50, # change this to how many ever iterations
    seed = round(runif(1) * 2^30) + 1:n()
  )



# Just checking!! ---------------------------------------------------------

quick_params <- params %>% 
  filter(batch %in% 1:3)

rm(design_factors, params)
source_obj <- ls()



#--------------------------------------------------------
# run simulations in parallel - mdply workflow
#--------------------------------------------------------

library(Pusto)

cluster <- start_parallel(source_obj = source_obj, 
                          setup = "register")



system.time(results <- plyr::mdply(quick_params, 
                                   .fun = run_sim, 
                                   .parallel = TRUE))

stop_parallel(cluster)



#--------------------------------------------------------
# Save results and details
#--------------------------------------------------------

session_info <- sessionInfo()
run_date <- date()

# batch names 
which_batches <- unique(quick_params$batch)
results_file <- paste0("sim_test_2_", paste(which_batches, collapse = "_"), ".RData")

# save
save(quick_params, results, session_info, run_date, file = results_file)