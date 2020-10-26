library(dplyr)
library(purrr)
library(mvtnorm)
library(robumeta)
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
                    test_dat,
                    R,
                    full_form = "X1 + X2 + X3 + X4 + X5", 
                    design_matrix = design_mat, 
                    seed = NULL) {

  
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
      full_model <- robu(as.formula(paste("g ~ ", full_form)), 
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
      cwb_params <- test_dat %>%
        select(null_model, indices_test) %>%
        mutate(R = R,
               full_form = full_form)
      
      boot_res <- pmap_dfr(cwb_params, 
                           .f = cwb, 
                           full_mod_org = full_model, 
                           cov_mat = cov_mat_cr1, 
                           dat = meta_data)
      
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
source_obj <- ls()

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

to_test_beta <- cross_df(list(beta_type = design_factors$beta_type, 
                              cov_test = to_test$cov_test)) %>%
  left_join(to_test %>% select(cov_test, null_model, indices_test), by = "cov_test") %>%
  mutate(keep = case_when(beta_type == "A" ~ TRUE,
                          beta_type %in% c("B1", "B5") ~ str_detect(cov_test, "X1"),
                          beta_type %in% c("C1", "C5") ~ str_detect(cov_test, "X2"),
                          beta_type %in% c("D1", "D5") ~ str_detect(cov_test, "X3"),
                          beta_type %in% c("E1", "E5") ~ str_detect(cov_test, "X4"),
                          beta_type %in% c("F1", "F5") ~ str_detect(cov_test, "X5"),
                          )) %>%
  filter(keep) %>%
  select(-keep) %>%
  group_by(beta_type) %>%
  nest() %>%
  ungroup()

params <- params %>%
  left_join(to_test_beta, by = "beta_type") %>%
  rename(test_dat = data)

glimpse(params)




# Just checking!! ---------------------------------------------------------

quick_params <- params %>% 
  filter(batch == 1) %>%
  mutate(R = 2,
         iterations = 2)

glimpse(quick_params)


system.time(
  results <-
    quick_params %>%
    mutate(res = pmap(., .f = run_sim)) %>%
    unnest(cols = res)
)

# 2248.625  user  27.416 system 2508.889 elapsed on 1023
# 1443.501  user  19.938 system 1528.106 elapsed on 1026

save(results, file = "../data/res_run_sim_1026.RData")


# FURRR -------------------------------------------------------------------


library(future)
library(furrr)
plan(multisession)

quick_params <- params %>% 
  filter(batch == 1) %>%
  mutate(R = 2,
         iterations = 2)

glimpse(quick_params)


# user 13.239   system 0.725 elapsed 795.245 

system.time(
  results_furrr <-
    quick_params %>%
    mutate(res = future_pmap(., .f = run_sim)) %>%
    unnest(cols = res)
)

save(results, file = "../data/res_run_sim_1026_parallel.RData")



# In parallel with furrr -------------------------------------------------

library(future)
library(furrr)
plan(multiprocess)

quick_params <- params %>% 
  filter(batch == 1) %>%
  mutate(R = 8,
         iterations = 4)

glimpse(quick_params)

system.time(
  results <-
    quick_params %>%
    mutate(res = future_pmap(., .f = run_sim)) %>%
    unnest(cols = res)
)

save(results, file = "../data/res_run_sim_1023_parallel.RData")



#--------------------------------------------------------
# run simulations in parallel - mdply workflow
#--------------------------------------------------------

library(Pusto)

cluster <- start_parallel(source_obj = source_obj)

system.time(results <- plyr::mdply(params, .fun = run_sim, .parallel = TRUE))

stopCluster(cluster)




#--------------------------------------------------------
# Save results and details
#--------------------------------------------------------

session_info <- sessionInfo()
run_date <- date()

save(params, results, session_info, run_date, file = "simulation_results_study_1.Rdata")