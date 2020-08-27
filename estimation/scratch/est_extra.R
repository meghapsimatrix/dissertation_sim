# extract residuals -------------------------------------------------------


change_to_mat <- function(res){
  
  as.matrix(res)
  
}


# extract the residuals

extract_res <- function(mod, dat = meta_data) {
  
  res <- clubSandwich:::residuals_CS.robu(mod)
  study <- dat$study
  
  split_res <- split(res, study)
  
  e_tilde_j <- map(split_res, change_to_mat)
  
  return(list(e_tilde_j))
  
}

system.time(null_res <- map(null_mods, extract_res))
null_res <- flatten(null_res)


# extract null cr models --------------------------------------------------

extract_B <- function(mod) {
  
  attr(vcovCR(mod, type = "CR2"), "adjustments")
  
}

extract_B(null_mods[[31]])  # works till 30
# 31 doesn't have a vcov adj matrix 
# because the model only has an intercept

null_mod <- robu(g ~ 1, 
                 studynum = study, 
                 var.eff.size = var_g,
                 small = FALSE,
                 data = meta_data)


attr(vcovCR(null_mod, type = "CR2"), "adjustments")

system.time(null_B <- map(null_mods[1:30], extract_B))

null_B

# multiply residuals by adj matrices --------------------------------------

#multiply B * e_tilde_j

mult_mat <- function(x, y) {
  
  as.vector(x %*% y)
  
}


t_res <- map2(.x = null_B, .y = null_res_30,
              ~ map2(.x, .y, mult_mat))


# res and t_res as vectors instead of matrices ----------------------------

null_cr2_res <- map(t_res, unlist)  # the names are all messed up?
null_res <- map(null_res_30, unlist)

head(transformed_res, 3)

# need to write out parts needed to run the cwb
# original_res
# transformed_res 
# type of test single or mch
# constraints - indices to test
# type of test 

# have all of these in a tibble or something then run cwb for each?

params <- tibble(indices = to_test$indices[-31],
                 null_res = null_res,
                 null_cr2_res = null_cr2_res)


# Extract stats for cwb ---------------------------------------------------

extract_stats <- function(mod, C, vcov_mat, method){
  
  Wald_test(mod, constraints = C, vcov = vcov_mat, test = "Naive-F") %>%
    as_tibble() %>%
    mutate(type = method)
}



# cluster wild bootstrapping ----------------------------------------------

cwb <- function(dat, indices){
  
  num_studies <- unique(dat$study)
  wts <- sample(c(-1, 1), size = length(num_studies), replace = TRUE)
  k_j <- as.numeric(table(dat$study))
  
  dat$eta <- rep(wts, k_j)
  
  dat$new_t <- with(dat, pred_null + res_null * eta)
  dat$new_t_adj <- with(dat, pred_null + t_res * eta)
  
  
  full_mod <- fit_mod("new_t ~ X1 + X2 + X3 + X4 + X5", dat)
  full_mod_adj <- fit_mod("new_t_adj ~ X1 + X2 + X3 + X4 + X5", dat)
  
  cov_mat <- vcovCR(full_mod, type = "CR1")
  cov_mat_adj <- vcovCR(full_mod_adj, type = "CR1")
  
  res <- extract_stats(full_mod, constrain_zero(indices), cov_mat, "CWB")
  res_adj <- extract_stats(full_mod_adj, constrain_zero(indices), cov_mat_adj, "CWB Adjusted")
  
  
  res <- bind_rows(res, res_adj)
  
  return(res)
  
}
