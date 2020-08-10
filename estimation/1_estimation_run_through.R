library(tidyverse)
library(robumeta)
library(clubSandwich)

# files -------------------------------------------------------------------
load("data/meta_data_practice.RData")

# combinations ------------------------------------------------------------

covs <- c("X1", "X2", "X3", "X4", "X5")

# terms combinations
comb_terms <- function(m, terms = covs) combn(length(terms), m, simplify = FALSE)

indices <- unlist(map(seq_along(1:5), comb_terms), recursive = FALSE)
equations <- map_chr(indices, function(x) paste(terms[x], collapse = "+"))

full_mod_indices <- 1:5


# null model 
null_indices <- map(indices, function(x) full_mod_indices[-x])
null_terms <- map(null_indices, function(x) terms[x])
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
full_formula <- "smd ~ X1 + X2 + X3 + X4 + X5"


full_mod <- robu(as.formula(full_formula), 
                 studynum = study, 
                 var.eff.size = var_smd,
                 small = FALSE,
                 data = dat)

full_mod

# function to run wald test -----------------------------------------------


get_wald <- function(model, constraints, cov_mat, test){
         
  Wald_test(model, 
            constraints = constrain_zero(constraints), 
            vcov = cov_mat,
            test = test)
}


# run wald ----------------------------------------------------------------

cov_mat_cr1 <- vcovCR(full_mod, type = "CR1")
cov_mat_cr2 <- vcovCR(full_mod, type = "CR2")

full_mod <- list(full_mod)
cov_mat_cr1 <- list(cov_mat_cr1)
cov_mat_cr2 <- list(cov_mat_cr2)
naive <- list("Naive-F")
htz <- list("HTZ")



system.time(naive_test <- pmap_dfr(list(full_mod, to_test$indices, cov_mat_cr1, "Naive-F"), get_wald))
system.time(htz_test <- pmap_dfr(list(full_mod, to_test$indices, cov_mat_cr2, "HTZ"), get_wald))

naive_res <- bind_cols(to_test %>% select(equ, type), naive_test)
htz_res <- bind_cols(to_test %>% select(equ, type), htz_test)

# CWB ---------------------------------------------------------------------

# fit null mods

fit_mod <- function(equation, dat = meta_data){
  
  robu(as.formula(equation), 
       studynum = study, 
       var.eff.size = var_smd,
       small = FALSE,
       data = dat)

}

params <- tibble(equation = paste("smd ~ ", to_test$null_model)) 
  
system.time(null_mods <- map(params$equation, fit_mod))


# extract residuals -------------------------------------------------------


change_to_mat <- function(res){
  
  as.matrix(res)
  
}




# extract the residuals

extract_res <- function(mod, dat = meta_data){

  res <- clubSandwich:::residuals_CS.robu(mod)
  study <- dat$study

  split_res <- split(res, study)

  e_tilde_j <- map(split_res, change_to_mat)

  return(list(e_tilde_j))

}

system.time(null_res <- map(null_mods, extract_res))
null_res_30 <- null_res[1:30]


# extract null cr models --------------------------------------------------

# this is not running

extract_B <- function(mod){
  
  attr(vcovCR(mod, type = "CR2"), "adjustments")
  
}

extract_B(null_mods[[31]])  # works till 30
# 31 doesn't have a vcov adj matrix 
# because the model only has an intercept

null_mod <- robu(smd ~ 1, 
                 studynum = study, 
                 var.eff.size = var_smd,
                 small = FALSE,
                 data = meta_data)


attr(vcovCR(null_mod, type = "CR2"), "adjustments")

system.time(null_B <- map(null_mods[1:30], extract_B))



# multiply residuals by adj matrices --------------------------------------

#multiply B * e_tilde_j

mult_mat <- function(x, y){
  
  x %*% y
  
}

# not working right now 
t_res <- map2(.x = null_B, .y = null_res_30,
              ~ map2(.x, .y, mult_mat))


# res and t_res as vectors instead of matrices ----------------------------



