library(tidyverse)
library(robumeta)
library(clubSandwich)

# files -------------------------------------------------------------------
load("data/meta_data_practice.RData")

# combinations ------------------------------------------------------------

covs <- c("X1", "X2", "X3", "X4", "X5", "X6_B+X6_C")

# terms combinations
comb_terms <- function(m, terms = covs) combn(length(terms), m, simplify = FALSE)

indices <- flatten(map(seq_along(1:6), comb_terms))
indices_test <- map(indices, function(x) x + 1)


indices_test_2 <- map(indices_test, function(x) ifelse(x == 7, c(7, 8), x))


equations <- map_chr(indices, function(x) paste(covs[x], collapse = "+"))

full_mod_indices <- 1:6

terms <- c("X1", "X2", "X3", "X4", "X5", "X6_B", "X6_C")

# null model 
null_indices <- map(indices, function(x) full_mod_indices[-x])
null_terms <- map(null_indices, function(x) covs[x])
null_model <- map_chr(null_terms, function(x) paste(x, collapse = "+"))


# all the null model to fit
to_test <- tibble(equ = equations,
                  #type = c(rep("single", 7), rep("mch", 26)),
                  indices = indices) %>%
  mutate(null_indices = null_indices,
         null_terms = null_terms, 
         null_model = null_model, 
         null_model = if_else(null_model == "", "1", null_model),
         indices_test = indices_test)


to_test <- to_test %>%
  mutate(indices_test_2 = ifelse(str_detect(equ, "X6"), c(indices_test, 8), indices_test))

# full model --------------------------------------------------------------

dat <- meta_data
full_formula <- "g ~ X1 + X2 + X3 + X4 + X5"


full_model <- robu(g ~ 0 + X1 + X2 + X3 + X4 + X5, 
                 studynum = study, 
                 var.eff.size = var_g,
                 small = FALSE,
                 data = dat)

full_model


# function to run wald test -----------------------------------------------

#a) fitted full model and b) null model specification/indices, and t
# he outputs would be a data frame with the test results corresponding to that hypothesis. 

estimate_wald <- function(model, indices_test, cov_mat, test){
  
  res <- Wald_test(model, 
                   constraints = constrain_zero(indices_test), 
                   vcov = cov_mat,
                   test = test) %>%
    select(test, p_val)
  
  return(res)
}


estimate_wald(full_model, indices_test = 2, cov_mat = vcovCR(full_model, type = "CR1"), test = "Naive-F")
estimate_wald(full_model, indices_test = 2, cov_mat = vcovCR(full_model, type = "CR1"), test = "HTZ")



# run wald ----------------------------------------------------------------

cov_mat_cr1 <- vcovCR(full_model, type = "CR1")
cov_mat_cr2 <- vcovCR(full_model, type = "CR2")

system.time(naive_test <- map_dfr(to_test$indices_test, estimate_wald, model = full_model, cov_mat = cov_mat_cr1, test = "Naive-F"))
system.time(htz_test <- map_dfr(to_test$indices_test, estimate_wald, model = full_model, cov_mat = cov_mat_cr2, test = "HTZ"))

naive_res <- bind_cols(to_test %>% select(equ, type), naive_test)
htz_res <- bind_cols(to_test %>% select(equ, type), htz_test)

# CWB ---------------------------------------------------------------------

# fit null mods

fit_mod <- function(equation, dat = meta_data) {
  
  robu(as.formula(equation), 
       studynum = study, 
       var.eff.size = var_g,
       small = FALSE,
       data = dat)

}

params <- tibble(equation = paste("g ~ ", to_test$null_model)) 
system.time(null_mods <- map(params$equation, fit_mod))




# run the cwb -------------------------------------------------------------

change_to_mat <- function(res){
  
  as.matrix(res)
  
}

mult_mat <- function(x, y) {
  
  as.vector(x %*% y)
  
}


extract_stats <- function(mod, C, vcov_mat, method){
  
  Wald_test(mod, constraints = C, vcov = vcov_mat, test = "Naive-F") %>%
    as_tibble() %>%
    mutate(test = method) %>%
    select(test, Fstat)
}



cwb <- function(dat, null_mod, full_mod = full_model, indices) {
  

  # residuals and transformed residuals -------------------------------------

  dat$res <- clubSandwich:::residuals_CS.robu(null_mod)
  dat$pred <- with(dat, g - res)
  split_res <- split(dat$res, dat$study)
  e_tilde_j <- map(split_res, change_to_mat)
  B_j <- attr(vcovCR(null_mod, type = "CR2"), "adjustments")
  dat$t_res <- unlist(pmap(list(B_j, e_tilde_j), mult_mat))
  

  # Rademacher weights ------------------------------------------------------
  
  num_studies <- unique(dat$study)
  k_j <- as.numeric(table(dat$study))
  
  set.seed(8062020)
  
  system.time(
    
  bootstraps <- rerun(.n = 399, {

  wts <- sample(c(-1, 1), size = length(num_studies), replace = TRUE)
  dat$eta <- rep(wts, k_j)
  dat$new_t <- with(dat, pred + res * eta)
  dat$new_t_adj <- with(dat, pred + t_res * eta)
  
  full_mod <- fit_mod("new_t ~ X1 + X2 + X3 + X4 + X5", dat)
  full_mod_adj <- fit_mod("new_t_adj ~ X1 + X2 + X3 + X4 + X5", dat)
  
  cov_mat <- vcovCR(full_mod, type = "CR1")
  cov_mat_adj <- vcovCR(full_mod_adj, type = "CR1")
  
  res <- extract_stats(full_mod, constrain_zero(indices), cov_mat, "CWB")
  res_adj <- extract_stats(full_mod_adj, constrain_zero(indices), cov_mat_adj, "CWB Adjusted")
  
  bind_rows(res, res_adj)
  
  }) %>%
    bind_rows()
  
  )
  
  org_F <- Wald_test(full_mod, 
                     constraints = constrain_zero(indices), 
                     vcov = vcovCR(full_mod, type = "CR1"),
                     test = "Naive-F") %>%
    pull(Fstat)
  
  
  p_boot <- bootstraps %>%
    group_by(test) %>%
    summarize(p_val = mean(Fstat > org_F)) %>%
    ungroup() 
  
  
  return(p_boot)
  
}


# Similar to HTZ and Naive F

system.time(res <- cwb(dat = dat, null_mod = null_mods[[1]], indices = 2))

Wald_test(full_model, constraints = constrain_zero(2), 
          vcov = vcovCR(full_model, type = "CR2"), test = "HTZ")

Wald_test(full_model, constraints = constrain_zero(2), 
          vcov = vcovCR(full_model, type = "CR1"), test = "Naive-F")
