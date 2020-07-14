library(tidyverse)
library(robumeta)
library(clubSandwich)
library(janitor)


# Load data and clean -----------------------------------------------------

load("data/tsl_dat.RData")

# complete cases 
tsl_dat_small <- tsl_dat %>%
  select(studyid, esid, delta, v, dv, g2age) %>%
  drop_na() %>%
  group_by(studyid) %>%
  mutate(k_j = n_distinct(esid),
         sigma_sq_j = mean(v),
         w_j = 1/(k_j * sigma_sq_j))


mod_prelim <- lm(delta ~ dv + g2age, data = tsl_dat_small, weights = w_j)
qe <- deviance(mod_prelim)



# calculate v -------------------------------------------------------------

full_x <- model.matrix(mod_prelim)

Xj <- full_x %>%
  as_tibble() %>%
  mutate(studyid = tsl_dat_small$studyid) %>%
  group_by(studyid) %>%
  clean_names()

X_j_mat <- group_split(Xj)

w_dat <- tsl_dat_small %>%
  select(studyid, w_j) %>%
  distinct(.)
  group_by(studyid)

wj <- group_split(w_dat)


change_to_mat <- function(dat){
  
  dat %>%
    select(-studyid) %>%
    as.matrix()
  
}

change_to_vec <- function(dat){
  
  dat %>%
    select(-studyid) %>%
    pull(w_j)
}


X_j <- map(X_j_mat, change_to_mat)
w_j <- map(wj, change_to_vec)

calc_v <- function(x, w){
  
  x_m <- t(x) %*% x

  w * x_m
  
} 


V_all <- pmap(list(X_j, w_j), calc_v)

V <- solve(V_all %>% reduce(`+`)) 

# wts for corr effects model
# wij = 1/ (kj * (sigma sq j + tau sq))

