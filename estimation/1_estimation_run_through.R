library(tidyverse)
library(robumeta)
library(clubSandwich)

# files -------------------------------------------------------------------
load("data/meta_data_practice.RData")

# combinations ------------------------------------------------------------

covs <- c("X1", "X2", "X3", "X4", "X5")

# terms combinations
comb_terms <- function(m, terms = covs) combn(length(terms), m, simplify = FALSE)

indices <- flatten(map(seq_along(1:5), comb_terms))
equations <- map_chr(indices, function(x) paste(covs[x], collapse = "+"))

full_mod_indices <- 1:5


# null model 
null_indices <- map(indices, function(x) full_mod_indices[-x])
null_terms <- map(null_indices, function(x) covs[x])
null_model <- map_chr(null_terms, function(x) paste(x, collapse = "+"))


# all the null model to fit
to_test <- tibble(equ = equations,
                  type = c(rep("single", 5), rep("mch", 26)),
                  indices = indices) %>%
  mutate(null_indices = null_indices,
         null_terms = null_terms, 
         null_model = null_model, 
         null_model = if_else(null_model == "", "1", null_model))


# full model --------------------------------------------------------------

dat <- meta_data
full_formula <- "g ~ X1 + X2 + X3 + X4 + X5"


full_mod <- robu(as.formula(full_formula), 
                 studynum = study, 
                 var.eff.size = var_g,
                 small = FALSE,
                 data = dat)

full_mod


# function to run wald test -----------------------------------------------

#a) fitted full model and b) null model specification/indices, and t
# he outputs would be a data frame with the test results corresponding to that hypothesis. 

get_wald <- function(model, indices, cov_mat, test){
         
  Wald_test(model, 
            constraints = constrain_zero(indices), 
            vcov = cov_mat,
            test = test)
}

get_wald(full_mod, indices = 2, cov_mat = vcovCR(full_mod, type = "CR1"), test = "Naive-F")
get_wald(full_mod, indices = 3, cov_mat = vcovCR(full_mod, type = "CR1"), test = "HTZ")


# run wald ----------------------------------------------------------------

cov_mat_cr1 <- vcovCR(full_mod, type = "CR1")
cov_mat_cr2 <- vcovCR(full_mod, type = "CR2")

system.time(naive_test <- map_dfr(to_test$indices, get_wald, model = full_mod, cov_mat = cov_mat_cr1, test = "Naive-F"))
system.time(htz_test <- map_dfr(to_test$indices, get_wald, model = full_mod, cov_mat = cov_mat_cr2, test = "HTZ"))

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
    mutate(type = method)
}



cwb <- function(dat, null_mod, indices) {
  

  # residuals and transformed residuals -------------------------------------

  dat$res <- clubSandwich:::residuals_CS.robu(null_mod)
  dat$pred <- with(dat, g - res)
  split_res <- split(dat$res, dat$study)
  e_tilde_j <- map(split_res, change_to_mat)
  B_j <- attr(vcovCR(null_mod, type = "CR2"), "adjustments")
  dat$t_res <- unlist(pmap(list(B_j, e_tilde_j), mult_mat))
  

  # Rademacher weights ------------------------------------------------------
  
  num_studies <- unique(dat$study)
  wts <- sample(c(-1, 1), size = length(num_studies), replace = TRUE)
  k_j <- as.numeric(table(dat$study))
  dat$eta <- rep(wts, k_j)
  dat$new_t <- with(dat, pred + res * eta)
  dat$new_t_adj <- with(dat, pred + t_res * eta)
  
  

  # fit the models ----------------------------------------------------------

  
  full_mod <- fit_mod("new_t ~ X1 + X2 + X3 + X4 + X5", dat)
  full_mod_adj <- fit_mod("new_t_adj ~ X1 + X2 + X3 + X4 + X5", dat)
  

  # extract test stats -------------------------------------------------------

  cov_mat <- vcovCR(full_mod, type = "CR1")
  cov_mat_adj <- vcovCR(full_mod_adj, type = "CR1")
  
  res <- extract_stats(full_mod, constrain_zero(indices), cov_mat, "CWB")
  res_adj <- extract_stats(full_mod_adj, constrain_zero(indices), cov_mat_adj, "CWB Adjusted")
  
  
  res <- bind_rows(res, res_adj)
  
  return(res)
  
}


cwb(dat = dat, null_mod = null_mods[[1]], indices = 2)

