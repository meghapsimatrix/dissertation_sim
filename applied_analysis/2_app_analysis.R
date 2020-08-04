library(tidyverse)
library(broom)
library(robumeta)
library(janitor)
library(xtable)
library(clubSandwich)

load("data/tsl_dat_20.RData")


robu_comp_tsl <- robu(delta ~ g2age + dv, 
                      studynum = study, 
                      var.eff.size = v,
                      small = TRUE,
                      data = tsl_dat)


# Single coefficient ------------------------------------------------------

res_age_naive <- Wald_test(robu_comp_tsl, 
          constraints = constrain_zero(2), 
          vcov = "CR0",
          test = "Naive-F")

res_age_htz <- Wald_test(robu_comp_tsl, 
          constraints = constrain_zero(2), 
          vcov = "CR2",
          test = "HTZ")




# multiple contrast hypothesis --------------------------------------------

robu_comp_tsl$reg_table$b.r
constrain_zero(3:7, robu_comp_tsl$reg_table$b.r)

res_mch_naive <- Wald_test(robu_comp_tsl, 
                            constraints = constrain_zero(3:7), 
                            vcov = "CR0", 
                            test = "Naive-F")

res_mch_htz <- Wald_test(robu_comp_tsl, 
                            constraints = constrain_zero(3:7), 
                            vcov = "CR2",
                            test = "HTZ")

res_mch <- bind_rows(res_mch_naive, res_mch_htz) %>%
  mutate(Method = rownames(.)) %>%
  select(Method, `F` = Fstat, delta, df_num, df_denom, p = p_val) %>%
  mutate_if(is.numeric, round, 3)

print(xtable(res_mch, digits = 3), file = "tsl_mch_app.tex")

      

# cluster wild ------------------------------------------------------------


full_mod <- robu(delta ~ g2age + dv, 
                         studynum = study, 
                         var.eff.size = v,
                         small = TRUE,
                         data = tsl_dat)

null_mod <- robu(delta ~ g2age, 
                      studynum = study, 
                      var.eff.size = v,
                      small = TRUE,
                      data = tsl_dat)

Wald_test(full_mod, vcov = "CR2", constraints = constrain_zero(2), test = "Naive-F")
Wald_test(full_mod, vcov = "CR2", constraints = constrain_zero(3:7), test = "Naive-F")

single_F <- Wald_test(full_mod, vcov = "CR2", test = "Naive-F", constraints = constrain_zero(2)) %>%
  as_tibble() %>%
  pull(Fstat)

mch_F <- Wald_test(full_mod, vcov = "CR2", test = "Naive-F", constraints = constrain_zero(3:7)) %>%
  as_tibble() %>%
  pull(Fstat)

tsl_dat <- tsl_dat %>%
  mutate(res_null = clubSandwich:::residuals_CS.robu(null_mod),
         pred_null = delta - res_null)


# adjustment matrix -------------------------------------------------------

B_j <- attr(vcovCR(full_mod, type = "CR2"), "adjustments")

e_tilde_j_all <- tsl_dat[, c("study", "res_null")]
e_tilde_j <- split(e_tilde_j_all, e_tilde_j_all$study)

change_to_mat <- function(dat){
  
  as.matrix(dat[, 2])
  
}

mult_mat <- function(x, y){
  
  x %*% y
  
}

e_tilde_j <- map(e_tilde_j, change_to_mat)


# transformed new residuals
tsl_dat$t_res <- unlist(pmap(list(B_j, e_tilde_j), mult_mat))


extract_stats <- function(mod, constraints, vcov_mat, method, var){
  
  Wald_test(mod, constraints = constraints, vcov = vcov_mat, test = "Naive-F") %>%
  as_tibble() %>%
  mutate(type = method,
         constraint = var)
}

# cov_mat <- vcovCR(full_mod, type = "CR2")
# extract_stats(full_mod, constrain_zero(2), cov_mat, "CWB", "age")

#extract_stats(robu_comp_tsl, constrain_zero(2), "CWB", "age")

# cluster wild bootstrapping ----------------------------------------------

cwb <- function(dat){
  
  num_studies <- unique(dat$study)
  wts <- sample(c(-1, 1), size = length(num_studies), replace = TRUE)
  k_j <- as.numeric(table(dat$study))
  
  dat$eta <- rep(wts, k_j)
    
  dat$new_t <- with(dat, pred_null + res_null * eta)
  dat$new_t_adj <- with(dat, pred_null + t_res * eta)
  
  full_mod <- robu(new_t ~ g2age + dv, 
                   studynum = study, 
                   var.eff.size = v,
                   small = FALSE,
                   data = dat)
  
  full_mod_adj <- robu(new_t_adj ~ g2age + dv, 
                       studynum = study, 
                       var.eff.size = v,
                       small = FALSE,
                       data = dat)
  
  cov_mat <- vcovCR(full_mod, type = "CR2")
  cov_mat_adj <- vcovCR(full_mod_adj, type = "CR2")
  
  res_single <- extract_stats(full_mod, constrain_zero(2), cov_mat, "CWB", "age")
  res_single_adj <- extract_stats(full_mod_adj, constrain_zero(2), cov_mat_adj, "CWB Adjusted", "age")
  res_mch <- extract_stats(full_mod, constrain_zero(3:7), cov_mat, "CWB", "mch")
  res_mch_adj <- extract_stats(full_mod_adj, constrain_zero(3:7), cov_mat_adj, "CWB Adjusted", "mch")
  
  res <- bind_rows(res_single, res_single_adj, res_mch, res_mch_adj)
  
  return(res)
  
}


cwb(tsl_dat)


set.seed(7232020)

system.time(

bootstraps <- rerun(.n = 999, {
  
  cwb(tsl_dat)
  
}) %>%
  bind_rows()

)


save(bootstraps, file = "data/bootstrap_reps.RData") 


calc_p_value <- function(dat = bootstraps, var, method, comp){
  
  dat %>%
    filter(constraint == var & type == method) %>%
    summarize(p = mean(Fstat > comp)) %>%
    mutate(Method = method)
  
}



res_age_cwb <- calc_p_value(var = "age", method = "CWB", comp = single_F)
res_age_cwb_adj <- calc_p_value(var = "age", method = "CWB Adjusted", comp = single_F)
res_mch_cwb <- calc_p_value(var = "mch", method = "CWB", comp = mch_F)
res_mch_cwb_adj <- calc_p_value(var = "mch", method = "CWB Adjusted", comp = mch_F)


# output the results ------------------------------------------------------

# singe coefs

res_sc <- bind_rows(res_age_naive, res_age_htz) %>%
  mutate(Method = test) %>%
  select(Method, `F` = Fstat, delta, df_num, df_denom, p = p_val) %>%
  bind_rows(res_age_cwb, res_age_cwb_adj) %>%
  mutate_if(is.numeric, round, 3)


print(xtable(res_sc, digits = 3), file = "tsl_sc_app.tex")


# multi-contrast hypothesis

res_mch <- bind_rows(res_mch_naive, res_mch_htz) %>%
  mutate(Method = test) %>%
  select(Method, `F` = Fstat, delta, df_num, df_denom, p = p_val) %>%
  bind_rows(res_mch_cwb, res_mch_cwb_adj) %>%
  mutate_if(is.numeric, round, 3)

print(xtable(res_mch, digits = 3), file = "tsl_mch_app.tex")