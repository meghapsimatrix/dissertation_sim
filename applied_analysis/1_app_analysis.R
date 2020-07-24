library(tidyverse)
library(broom)
library(robumeta)
library(janitor)
library(xtable)
library(clubSandwich)

load("data/tsl_dat.RData")

tsl_dat <- tsl_dat %>%
  select(study, es_num, delta, v, dv, g2age) %>%
  drop_na() %>%
  mutate(dv = factor(dv),
         dv = fct_recode(dv, bac = "Blood alcohol concentration",
                             comb = "Combined measures (e.g., AUDIT)",
                             fhu = "Frequency of heavy use",
                             fu = "Frequency of use",
                             pc = "Peak consumption",
                             qu = "Quantity of use"))
  


robu_comp_tsl <- robu(delta ~ g2age + dv, 
                      studynum = study, 
                      var.eff.size = v,
                      small = TRUE,
                      data = tsl_dat)


# Single coefficient ------------------------------------------------------

res_age_naive <- Wald_test(robu_comp_tsl, 
          constraints = 2, 
          vcov = "CR0",
          test = "Naive-F")

res_age_htz <- Wald_test(robu_comp_tsl, 
          constraints = 2, 
          vcov = "CR2",
          test = "HTZ")




# multiple contrast hypothesis --------------------------------------------



res_mch_naive <- Wald_test(robu_comp_tsl, 
                            constraints = 3:7, 
                            vcov = "CR0", 
                            test = "Naive-F")

res_mch_htz <- Wald_test(robu_comp_tsl, 
                            constraints = 3:7, 
                            vcov = "CR2",
                            test = "HTZ")

res_mch <- bind_rows(res_mch_naive, res_mch_htz) %>%
  mutate(Method = rownames(.)) %>%
  select(Method, `F` = Fstat, delta, df, p = p_val) %>%
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

Wald_test(full_mod, vcov = "CR2", constraints = 2, test = "chi-sq")
Wald_test(full_mod, vcov = "CR2", constraints = 3:7, test = "chi-sq")

single_F <- Wald_test(full_mod, vcov = "CR2", test = "chi-sq", constraints = 2) %>%
  as_tibble() %>%
  pull(Fstat)

mch_F <- Wald_test(full_mod, vcov = "CR2", test = "chi-sq", constraints = 3:7) %>%
  as_tibble() %>%
  pull(Fstat)

tsl_dat <- tsl_dat %>%
  mutate(res_null = clubSandwich:::residuals_CS.robu(null_mod),
         pred_null = delta - res_null)


# adjustment matrix -------------------------------------------------------

B_j <- attr(vcovCR(null_mod, type = "CR2"), "adjustments")

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


extract_stats <- function(mod, constraints, method, var){
  
  Wald_test(mod, constraints = constraints, vcov = "CR2", test = "Naive-F") %>%
  as_tibble() %>%
  mutate(type = method,
         constraint = var)
}

extract_stats(robu_comp_tsl, 2, "CWB", "age")

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
                   small = TRUE,
                   data = dat)
  
  full_mod_adj <- robu(new_t_adj ~ g2age + dv, 
                       studynum = study, 
                       var.eff.size = v,
                       small = TRUE,
                       data = dat)
  
  res_single <- extract_stats(full_mod, 2, "CWB", "age")
  res_single_adj <- extract_stats(full_mod_adj, 2, "CWB Adjusted", "age")
  res_mch <- extract_stats(full_mod, 3:7, "CWB", "mch")
  res_mch_adj <- extract_stats(full_mod_adj, 3:7, "CWB Adjusted", "mch")
  
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

6777.061 /60

res_age_cwb <- calc_p_value(var = "age", method = "CWB", comp = single_F)
res_age_cwb_adj <- calc_p_value(var = "age", method = "CWB Adjusted", comp = single_F)
res_mch_cwb <- calc_p_value(var = "mch", method = "CWB", comp = mch_F)
res_mch_cwb_adj <- calc_p_value(var = "mch", method = "CWB Adjusted", comp = mch_F)


# output the results ------------------------------------------------------

# singe coefs

res_sc <- bind_rows(res_age_naive, res_age_htz) %>%
  mutate(Method = rownames(.)) %>%
  select(Method, `F` = Fstat, delta, df, p = p_val) %>%
  bind_rows(res_age_cwb, res_age_cwb_adj) %>%
  mutate_if(is.numeric, round, 3)


print(xtable(res_sc, digits = 3), file = "tsl_sc_app.tex")


# multi-contrast hypothesis

res_mch <- bind_rows(res_mch_naive, res_mch_htz) %>%
  mutate(Method = rownames(.)) %>%
  select(Method, `F` = Fstat, delta, df, p = p_val) %>%
  bind_rows(res_mch_cwb, res_mch_cwb_adj) %>%
  mutate_if(is.numeric, round, 3)

print(xtable(res_mch, digits = 3), file = "tsl_mch_app.tex")
