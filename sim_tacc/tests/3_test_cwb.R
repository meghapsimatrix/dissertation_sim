
null_model <- test_dat$null_model[[1]]
dat <- meta_data

null_mod <- robu(as.formula(null_model), 
                 studynum = study, 
                 var.eff.size = var_g,
                 small = FALSE,
                 data = dat)

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


    
    wts <- sample(c(-1, 1), size = length(num_studies), replace = TRUE)
    dat$eta <- rep(wts, k_j)
    dat$new_t <- with(dat, pred + res * eta)
    dat$new_t_adj <- with(dat, pred + t_res * eta)
    
    full_mod_form <- "X1 + X2 + X3 + X4 + X5"
    
    full_mod <- robu(as.formula(paste("new_t ~", full_mod_form)), 
                     studynum = study, 
                     var.eff.size = var_g,
                     small = FALSE,
                     data = dat)
    
    full_mod_adj <- robu(as.formula(paste("new_t_adj ~", full_mod_form)), 
                         studynum = study, 
                         var.eff.size = var_g,
                         small = FALSE,
                         data = dat)
    
    
    cov_mat <- vcovCR(full_mod, type = "CR1")
    cov_mat_adj <- vcovCR(full_mod_adj, type = "CR1")
    
    indices_test <- test_dat$indices_test[[1]]
    
    res <- extract_stats(full_mod, constrain_zero(indices_test), cov_mat, "CWB")
    res_adj <- extract_stats(full_mod_adj, constrain_zero(indices_test), cov_mat_adj, "CWB Adjusted")
    