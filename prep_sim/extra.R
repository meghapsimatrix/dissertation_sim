# Calculate CR0 -----------------------------------------------------------

dat <- meta_data
study <- dat$study
m <- length(unique(dat$study))
res <- update_fit$residuals


e_j <- split(res, study)
e_j <- map(e_j, change_to_mat)


X_j <- split(as.data.frame(X), study)
X_j <- map(X_j, change_to_mat)

k_j <- as.numeric(table(dat$study))


I_j <- map(k_j, diag)

w_j <- data.frame(study, w_j = update_fit$weights) %>%
  distinct(.) %>%
  pull(w_j)

W_j <- pmap(list(w_j, I_j), mult_mat_const)
W <- update_fit$weights * diag(nrow(dat))



M <- solve(t(X) %*% W %*% X)

meat <- pmap(list(X_j, W_j, e_j), multiply_rve)
meat_sum <- reduce(meat, `+`)
meat_sum <- as.matrix(meat_sum)

str(meat_sum)





# non zero betas only for power -------------------------------------------

if (beta_type %in% c("B1", "B5")){
  
  test_dat <- test_dat %>%
    filter(str_detect(cov_test, "X1"))
  
} else if(beta_type %in% c("C1", "C5")){
  
  test_dat <- test_dat %>%
    filter(str_detect(cov_test, "X2"))
  
} else if(beta_type %in% c("D1", "D5")){
  
  test_dat <- test_dat %>%
    filter(str_detect(cov_test, "X3"))
  
} else if(beta_type %in% c("E1", "E5")){
  
  test_dat <- test_dat %>%
    filter(str_detect(cov_test, "X4"))
  
} else if(beta_type %in% c("F1", "F5")){
  
  test_dat <- test_dat %>%
    filter(str_detect(cov_test, "X5"))
  
} else if (beta_type == "A" ){
  
  test_dat <- test_dat 
  
} 
