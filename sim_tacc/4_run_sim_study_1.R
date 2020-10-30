library(dplyr)
library(purrr)
library(mvtnorm)
library(clubSandwich)
library(tidyr)
library(stringr)

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

run_sim <- function(iterations, 
                    m, 
                    tau, 
                    rho, 
                    beta_type, 
                    batch,
                    R,
                    full_form = "X1 + X2 + X3 + X4 + X5", 
                    test_dat = to_test,
                    design_matrix = design_mat, 
                    seed = NULL) {
  
  require(dplyr)
  require(purrr)
  require(mvtnorm)
  require(clubSandwich)
  require(tidyr)
  require(stringr)
  
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
      cwb_params <- test_dat %>%
        select(null_model, indices_test)
      
      boot_res <- pmap_dfr(cwb_params, 
                           .f = cwb, 
                           dat = meta_data,
                           R = R,
                           full_form = full_form)
      
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

set.seed(20150316) # change this seed value!

# now express the simulation parameters as vectors/lists

design_factors <- list(
  m = c(10, 20, 40, 80), 
  tau = c(0.1, 0.3),
  rho = c(0.5, 0.8),
  R = 399,
  beta_type = c("A", "B1", "B5", "C1", "C5", "D1", "D5", "E1", "E5", "F1", "F5"),
  batch = 1:10
)

# combine into a design set


params <-
  cross_df(design_factors) %>%
  mutate(
    iterations = 100, # change this to how many ever iterations
    seed = round(runif(1) * 2^30) + 1:n()
  )

glimpse(params)

# to_test_beta <- cross_df(list(beta_type = design_factors$beta_type, 
#                               cov_test = to_test$cov_test)) %>%
#   left_join(to_test %>% select(cov_test, null_model, indices_test), by = "cov_test") %>%
#   mutate(keep = case_when(beta_type == "A" ~ TRUE,
#                           beta_type %in% c("B1", "B5") ~ str_detect(cov_test, "X1"),
#                           beta_type %in% c("C1", "C5") ~ str_detect(cov_test, "X2"),
#                           beta_type %in% c("D1", "D5") ~ str_detect(cov_test, "X3"),
#                           beta_type %in% c("E1", "E5") ~ str_detect(cov_test, "X4"),
#                           beta_type %in% c("F1", "F5") ~ str_detect(cov_test, "X5"),
#                           )) %>%
#   filter(keep) %>%
#   select(-keep) %>%
#   group_by(beta_type) %>%
#   nest() %>%
#   ungroup()
# 
# params <- params %>%
#   left_join(to_test_beta, by = "beta_type") %>%
#   rename(test_dat = data)
# 
# glimpse(params)
# 
# 


# Just checking!! ---------------------------------------------------------

quick_params <- params %>% 
  filter(batch == 1) %>%
  mutate(R = 2,
         iterations = 2)

glimpse(quick_params)

rm(design_factors, params)
source_obj <- ls()

source_obj

system.time(
  results <-
    quick_params %>%
    mutate(res = pmap(., .f = run_sim)) %>%
    unnest(cols = res)
)

# 2508.889 elapsed on 1023
# 1528.106 elapsed on 1026
# 1442.405 elapsed on 1026 afternoon
# 566.618 elapsed on 1027
# 323.27 elapsed on 1028 on R 3.6.0 on windows desktop
# 898.874 on 1030 on mac

save(results, file = "../data/res_run_sim_1030.RData")



# In parallel with furrr -------------------------------------------------

# library(future)
# library(furrr)
# 
# #plan(multiprocess)
# plan(multicore)
# 
# quick_params <- params %>% 
#   filter(batch == 1) %>%
#   mutate(R = 8,
#          iterations = 4)
# 
# glimpse(quick_params)
# 
# system.time(
#   results <-
#     quick_params %>%
#     mutate(res = future_pmap(., .f = run_sim)) %>%
#     unnest(cols = res)
# )
# 
# save(results, file = "../data/res_run_sim_1023_parallel.RData")


# Error in terms.default(object) : no terms component nor attribute
# Timing stopped at: 0.34 0.2 3.62
# furrr is not running for me in either windows

# Error: Problem with `mutate()` input `res`.
# x missing or negative weights not allowed
# ℹ Input `res` is `future_pmap(., .f = run_sim)`.

#--------------------------------------------------------
# run simulations in parallel - mdply workflow
#--------------------------------------------------------


library(Pusto)

cluster <- start_parallel(source_obj = source_obj, 
                          setup = "register",
                          packages = c("dplyr", "purrr", "mvtnorm",
                                       "clubSandwich",
                                       "tidyr", "stringr"))


# any function or object that R needs to have 

system.time(results <- plyr::mdply(quick_params, 
                                   .fun = run_sim, 
                                   .parallel = TRUE))

stopCluster(cluster)

# 446.282 on 1030

save(results, file = "../data/res_run_sim_mdply_1030.RData")

# Error in do.ply(i) : 
#   task 1 failed - "You must specify at least one constraint."
# Timing stopped at: 34.71 0.778 19.61

#--------------------------------------------------------
# Save results and details
#--------------------------------------------------------

session_info <- sessionInfo()
run_date <- date()

save(params, results, session_info, run_date, file = "simulation_results_study_1.Rdata")