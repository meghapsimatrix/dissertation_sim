
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
