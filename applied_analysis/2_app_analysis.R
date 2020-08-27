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
          vcov = "CR1",
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
                            vcov = "CR1", 
                            test = "Naive-F")

res_mch_htz <- Wald_test(robu_comp_tsl, 
                            constraints = constrain_zero(3:7), 
                            vcov = "CR2",
                            test = "HTZ")

Wald_test(robu_comp_tsl, 
          constraints = constrain_zero(3:7), 
          vcov = "CR2",
          test = "Naive-F")
      

# cluster wild ------------------------------------------------------------

full_mod <- robu(delta ~ g2age + dv, 
                 studynum = study, 
                 var.eff.size = v,
                 small = TRUE,
                 data = tsl_dat)

null_mod_mch <- robu(delta ~ g2age, 
                     studynum = study, 
                     var.eff.size = v,
                     small = TRUE,
                     data = tsl_dat)

null_mod_single <- robu(delta ~ dv, 
                        studynum = study, 
                        var.eff.size = v,
                        small = TRUE,
                        data = tsl_dat)



# save the F from the full model ------------------------------------------

Wald_test(full_mod, vcov = "CR1", constraints = constrain_zero(2), test = "Naive-F")
Wald_test(full_mod, vcov = "CR1", constraints = constrain_zero(3:7), test = "Naive-F")

single_F <- Wald_test(full_mod, vcov = "CR1", test = "Naive-F", constraints = constrain_zero(2)) %>%
  as_tibble() %>%
  pull(Fstat)

mch_F <- Wald_test(full_mod, vcov = "CR1", test = "Naive-F", constraints = constrain_zero(3:7)) %>%
  as_tibble() %>%
  pull(Fstat)



# Extract stats for cwb ---------------------------------------------------

extract_stats <- function(mod, constraints, vcov_mat, method){
  
  Wald_test(mod, constraints = constraints, vcov = vcov_mat, test = "Naive-F") %>%
    as_tibble() %>%
    mutate(type = method)
}



# cluster wild bootstrapping ----------------------------------------------

cwb <- function(dat, constraints) {
  
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
  
  cov_mat <- vcovCR(full_mod, type = "CR1")
  cov_mat_adj <- vcovCR(full_mod_adj, type = "CR1")
  
  res <- extract_stats(full_mod, constraints, cov_mat, "CWB")
  res_adj <- extract_stats(full_mod_adj, constraints, cov_mat_adj, "CWB Adjusted")
  
  res <- bind_rows(res, res_adj)
  
  return(res)
  
}




# functions for matrices --------------------------------------------------

change_to_mat <- function(dat){
  
  as.matrix(dat)
  
}

mult_mat <- function(x, y){
  
  x %*% y
  
}


# MCH ---------------------------------------------------------------------

B_j <- attr(vcovCR(null_mod_mch, type = "CR2"), "adjustments")


tsl_dat <- tsl_dat %>%
  mutate(res_null = clubSandwich:::residuals_CS.robu(null_mod_mch),
         pred_null = delta - res_null)

e_tilde_j <- split(tsl_dat$res_null, tsl_dat$study)
e_tilde_j <- map(e_tilde_j, change_to_mat)

# transformed new residuals
tsl_dat$t_res <- unlist(pmap(list(B_j, e_tilde_j), mult_mat))



cwb(tsl_dat, constraints = constrain_zero(2))


set.seed(7232020)

system.time(

  bootstraps_mch <- rerun(.n = 999, {
    
    cwb(tsl_dat, constraints = constrain_zero(3:7))
    
  }) %>%
    bind_rows()

)


save(bootstraps_mch, file = "data/bootstrap_mch.RData") 



# FOR SINGLE COEFS --------------------------------------------------------

B_j <- attr(vcovCR(null_mod_single, type = "CR2"), "adjustments")


tsl_dat <- tsl_dat %>%
  mutate(res_null = clubSandwich:::residuals_CS.robu(null_mod_single),
         pred_null = delta - res_null)

e_tilde_j <- split(tsl_dat$res_null, tsl_dat$study)
e_tilde_j <- map(e_tilde_j, change_to_mat)

# transformed new residuals
tsl_dat$t_res <- unlist(pmap(list(B_j, e_tilde_j), mult_mat))

set.seed(8062020)

system.time(
  
  bootstraps_age <- rerun(.n = 999, {
    
    cwb(tsl_dat, constraints = constrain_zero(2))
    
  }) %>%
    bind_rows()
  
)

save(bootstraps_age, file = "data/bootstrap_age.RData") 



# Calculate the bootstrap p value -----------------------------------------


age_boot <- bootstraps_age %>%
  group_by(type) %>%
  summarize(p = mean(Fstat > single_F)) %>%
  ungroup() %>%
  rename(Method = type)

mch_boot <- bootstraps_mch %>%
  group_by(type) %>%
  summarize(p = mean(Fstat > mch_F)) %>%
  ungroup() %>%
  rename(Method = type)


# output the results ------------------------------------------------------

# singe coefs

res_sc <- bind_rows(res_age_naive, res_age_htz) %>%
  select(Method = test, `F` = Fstat, delta, df_num, df_denom, p = p_val) %>%
  bind_rows(age_boot) %>%
  mutate_if(is.numeric, round, 3)


print(xtable(res_sc, digits = 3), file = "tsl_sc_app.tex")


# multi-contrast hypothesis

res_mch <- bind_rows(res_mch_naive, res_mch_htz) %>%
  select(Method = test, `F` = Fstat, delta, df_num, df_denom, p = p_val) %>%
  bind_rows(mch_boot) %>%
  mutate_if(is.numeric, round, 3)

print(xtable(res_mch, digits = 3), file = "tsl_mch_app.tex")
