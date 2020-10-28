library(robumeta)
library(clubSandwich)
library(microbenchmark)
library(tidyverse)

calculate_F <- function(beta, vcov, constraints, p = 6, test){
  
  C_mat <- diag(1L, nrow = p)[constraints,,drop = FALSE]    
  
  Q <- as.numeric(t(C_mat %*% beta) %*% chol2inv(chol(C_mat %*% vcov %*% t(C_mat))) %*% C_mat %*% beta)
  q <- nrow(C_mat)
  
  Fstat <- Q/q
  
  test_res <- tibble(test = test, 
                     Fstat = Fstat)
  
  return(test_res)
  
}


full_form <- "X1 + X2 + X3 + X4 + X5"
load("data/meta_data_practice.RData")


check <- robu(as.formula(paste("g ~ ", full_form)), 
              studynum = study, 
              var.eff.size = var_g,
              small = FALSE,
              data = meta_data)

check$reg_table$b.r
cov_mat_cr1 <- vcovCR(check, type = "CR1")
calculate_F(beta = check$reg_table$b.r, vcov = cov_mat_cr1, constraints = 3, test = "CWB Check")

Wald_test(check, 
          constraints = constrain_zero(3),
          vcov = cov_mat_cr1,
          test = "Naive-F", 
          tidy = TRUE) 



microbenchmark(
  calculate_F(beta = check$reg_table$b.r, vcov = cov_mat_cr1, constraints = 3:6, test = "CWB Check"),
  Wald_test(check, 
            constraints = constrain_zero(3:6),
            vcov = cov_mat_cr1,
            test = "Naive-F", 
            tidy = TRUE) 
)
