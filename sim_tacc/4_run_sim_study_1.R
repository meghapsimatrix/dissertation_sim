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
                    test_dat,
                    R = 399, # need to figure out how to pass this to cwb()
                    full_form = "X1 + X2 + X3 + X4 + X5", 
                    design_matrix = design_mat, 
                    seed = NULL) {

  
  if (!is.null(seed)) set.seed(seed)

  
  results <-
    map_dfr(1:iterations, {
        
      # generate data ------------------------------------------------------------
      meta_data <- generate_rmeta(m = m, 
                                  tau = tau, 
                                  rho = rho,
                                  covs = design_matrix,
                                  beta_type = beta_type)
      
      # Fit full model on data --------------------------------------------------
      full_formula <- full_form
      full_model <- robu(as.formula(paste("g ~ ", full_formula)), 
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
        select(null_model, indices_test)
      
      # if i don't put data and R and full_mod_form as default something goes wrong
      # need to figure out how to do R 
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
  R = 399,
  beta_type = c("A", "B1", "B5", "C1", "C5", "D1", "D5", "E1", "E5", "F1", "F5")
)

# combine into a design set

params <-
  cross_df(design_factors) %>%
  mutate(
    iterations = 1000, # change this to how many ever iterations
    seed = round(runif(1) * 2^30) + 1:n()
  )

to_test_beta <- cross_df(list(beta_type = design_factors$beta_type, 
                              cov_test = test_dat$cov_test)) %>%
  left_join(test_dat %>% select(cov_test, null_model, indices_test), by = "cov_test") %>%
  mutate(keep = case_when(beta_type == "A" ~ TRUE,
                          beta_type %in% c("B1", "B5") ~ str_detect(cov_test, "X1"),
                          beta_type %in% c("C1", "C5") ~ str_detect(cov_test, "X2"),
                          beta_type %in% c("D1", "D5") ~ str_detect(cov_test, "X3"),
                          beta_type %in% c("E1", "E5") ~  str_detect(cov_test, "X4"),
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