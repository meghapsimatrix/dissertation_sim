library(clubSandwich)
library(dplyr)
library(robumeta)

# Issues
# any way to extract the formula from robu?
# want to make this work for lm - what's the best way ? class
# what other types of models?
# output results as tibble or print or both
# Include the multiplication by CR2 or not
# better name for the function



cwb <- function(data,
                smd,
                var,
                cluster,
                full_model,
                null_model, 
                indices, 
                R = 999) {
  
  data$smd <- data %>%
    pull({{smd}})
  
  data$var <- data %>%
    pull({{var}})
  
  data$cluster <- data %>%
    pull({{cluster}})
  
  
  # residuals and transformed residuals -------------------------------------
  
  data$res <- clubSandwich:::residuals_CS.robu(null_model)
  data$pred <- with(data, smd - res)
  split_res <- split(data$res, data$cluster)

  
  # Rademacher weights ------------------------------------------------------
  
  num_cluster <- unique(data$cluster)
  k_j <- as.numeric(table(data$cluster))
  
  system.time(
    
    bootstraps <- rerun(.n = R, {
      
      wts <- sample(c(-1, 1), size = length(num_cluster), replace = TRUE)
      data$eta <- rep(wts, k_j)
      data$new_y <- with(data, pred + res * eta)

      boot_mod <- robu(as.formula(paste("new_y ~ ", 
                                  paste(full_model$reg_table$labels[!str_detect(full_model$reg_table$labels, "Intercept")], 
                                        collapse = "+"))), 
                       studynum = cluster, 
                       var.eff.size = var,
                       small = FALSE,
                       data = data)
      
      cov_mat <- vcovCR(boot_mod, type = "CR1")

      res <- Wald_test(boot_mod, 
                       constraints = constrain_zero(indices),
                       vcov = cov_mat, 
                       test = "Naive-F") 
      
    }) %>%
      bind_rows()
    
  )
  
  org_F <- Wald_test(full_model, 
                     constraints = constrain_zero(indices), 
                     vcov = vcovCR(full_model, type = "CR1"),
                     test = "Naive-F") %>%
    pull(Fstat)
  
  
  p_boot <- bootstraps %>%
    summarize(p_val = mean(Fstat > org_F)) %>%
    ungroup() %>%
    mutate(test = "CWB") %>%
    dplyr::select(test, p_val)
  
  
  return(p_boot)
  
}


load("data/meta_data_practice.RData")

full <- robu(g ~ X1 + X2 + X3 + X4 + X5, 
                   studynum = study, 
                   var.eff.size = var_g,
                   small = FALSE,
                   data = meta_data)

null <- robu(g ~ 1, 
             studynum = study, 
             var.eff.size = var_g,
             small = FALSE,
             data = meta_data)

# 
# data <- meta_data
# data$smd <- data$g
# data$var <- data$var_g
# data$cluster <- data$study
# null_model <- null
# full_model <- full
# 

system.time(res <- cwb(data = meta_data, 
    smd = g, 
    var = var_g, 
    cluster = study,
    full_model = full,
    null_model = null, 
    indices = 2:6,
    R = 999))
