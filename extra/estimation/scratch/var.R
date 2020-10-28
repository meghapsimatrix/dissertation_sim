library(tidyverse)
library(broom)
library(robumeta)
library(janitor)
library(xtable)
library(clubSandwich)

load("data/meta_data_practice.RData")

full_mod <- robu(g ~ X1 + X2 + X3 + X4 + X5, 
                 studynum = study, 
                 var.eff.size = var_g,
                 small = TRUE,
                 data = meta_data)

null_mod_mch <- robu(g ~ X1, 
                 studynum = study, 
                 var.eff.size = var_g,
                 small = TRUE,
                 data = meta_data)

# save the F from the full model ------------------------------------------

Wald_test(full_mod, vcov = "CR1", constraints = constrain_zero(3:6), test = "Naive-F")

mch_F <- Wald_test(full_mod, vcov = "CR1", test = "Naive-F", constraints = constrain_zero(3:6)) %>%
  as_tibble() %>%
  pull(Fstat)



# Extract stats for cwb ---------------------------------------------------

extract_stats <- function(mod, constraints, vcov_mat, method, var){
  
  Wald_test(mod, constraints = constraints, vcov = vcov_mat, test = "Naive-F") %>%
    as_tibble() %>%
    mutate(type = method,
           constraint = var)
}



# cluster wild bootstrapping ----------------------------------------------

cwb <- function(dat, single){
  
  num_studies <- unique(dat$study)
  wts <- sample(c(-1, 1), size = length(num_studies), replace = TRUE)
  k_j <- as.numeric(table(dat$study))
  
  dat$eta <- rep(wts, k_j)
  
  dat$new_t <- with(dat, pred_null + res_null * eta)
  dat$new_t_adj <- with(dat, pred_null + t_res * eta)
  
  
  full_mod <- robu(new_t ~ X1 + X2 + X3 + X4 + X5, 
                   studynum = study, 
                   var.eff.size = var_g,
                   small = FALSE,
                   data = dat)
  
  full_mod_adj <- robu(new_t_adj ~ X1 + X2 + X3 + X4 + X5, 
                       studynum = study, 
                       var.eff.size = var_g,
                       small = FALSE,
                       data = dat)
  
  cov_mat <- vcovCR(full_mod, type = "CR1")
  cov_mat_adj <- vcovCR(full_mod_adj, type = "CR1")
  
  if(single == TRUE){
    res <- extract_stats(full_mod, constrain_zero(2), cov_mat, "CWB", "age")
    res_adj <- extract_stats(full_mod_adj, constrain_zero(2), cov_mat_adj, "CWB Adjusted", "single")
  }
  else{
    res <- extract_stats(full_mod, constrain_zero(3:6), cov_mat, "CWB", "mch")
    res_adj <- extract_stats(full_mod_adj, constrain_zero(3:6), cov_mat_adj, "CWB Adjusted", "mch")
  }
  
  res <- bind_rows(res, res_adj)
  
  return(res)
  
}




# functions for matrices --------------------------------------------------

change_to_mat <- function(dat){
  
  as.matrix(dat[, 2])
  
}

mult_mat <- function(x, y){
  
  x %*% y
  
}


# MCH ---------------------------------------------------------------------

B_j <- attr(vcovCR(null_mod_mch, type = "CR2"), "adjustments")


meta_data <- meta_data %>%
  mutate(res_null = clubSandwich:::residuals_CS.robu(null_mod_mch),
         pred_null = g - res_null)


e_tilde_j_all <- meta_data[, c("study", "res_null")]
e_tilde_j <- split(e_tilde_j_all, e_tilde_j_all$study)
e_tilde_j <- map(e_tilde_j, change_to_mat)

# transformed new residuals
meta_data$t_res <- unlist(pmap(list(B_j, e_tilde_j), mult_mat))


cwb(meta_data, single = FALSE)


set.seed(7232020)

system.time(
  
  bootstraps_mch <- rerun(.n = 399, {
    
    cwb(meta_data, single = FALSE)
    
  }) %>%
    bind_rows()
  
)

# around 2 and a half mins with 999
# around 

dat <- meta_data
num_studies <- unique(dat$study)
wts <- sample(c(-1, 1), size = length(num_studies), replace = TRUE)
k_j <- as.numeric(table(dat$study))

dat$eta <- rep(wts, k_j)

dat$new_t <- with(dat, pred_null + res_null * eta)
dat$new_t_adj <- with(dat, pred_null + t_res * eta)

dat$d <- with(dat, g / (1 - (3/((4 * (N - 2)) - 1))))
dat$vd <- with(dat, (4 / N + d^2 / (2 * (N - 2))))
dat$vg <- with(dat, vd * (1 - (3/((4 * (N - 2)) - 1)))^2)

dat %>%
  select(vd, vg) %>%
  View()
