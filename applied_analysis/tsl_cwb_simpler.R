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
         w_tilde_j = 1/(k_j * sigma_sq_j)) %>%
  ungroup()


mod_prelim <- lm(delta ~ dv + g2age, data = tsl_dat_small, weights = w_tilde_j)


length(residuals(mod_prelim))

tsl_dat_small <- tsl_dat_small %>%
  mutate(res = residuals(mod_prelim))

QE <- with(tsl_dat_small, sum(w_tilde_j * res^2))

# calculate v -------------------------------------------------------------

full_x <- model.matrix(mod_prelim)

Xj <- full_x %>%
  as_tibble() %>%
  mutate(studyid = tsl_dat_small$studyid) %>%
  group_by(studyid) %>%
  clean_names()

X_j_mat <- group_split(Xj)

w_tilde_dat <- tsl_dat_small %>%
  select(studyid, w_tilde_j) %>%
  distinct(.) %>%
  group_by(studyid)

w_tilde_j <- group_split(w_tilde_dat)


num_es <- tsl_dat_small %>%
  select(studyid, k_j) %>%
  distinct(.) %>%
  select(k_j)

create_identity <- function(k){
  
  diag(k)
  
}


k <- num_es %>% pull(k_j)

I_j <- map(k, create_identity)

change_to_mat <- function(dat){
  
  dat %>%
    select(-studyid) %>%
    as.matrix()
  
}

change_to_vec <- function(dat){
  
  dat %>%
    select(-studyid) %>%
    pull(w_tilde_j)
}


X_j <- map(X_j_mat, change_to_mat)

w_tilde_j <- map(w_tilde_j, change_to_vec)

create_big_W <- function(w, id){
  
  w * id
  
}

W_tilde_j <- pmap(list(w_tilde_j, I_j), create_big_W)


create_m_tilde <- function(X, W){
  
  t(X) %*% W %*% X
  
}

M_tilde_all <- pmap(list(X_j, W_tilde_j), create_m_tilde)

M_tilde <- M_tilde_all %>% 
  reduce(`+`) %>%
  solve()


calculate_p <- function(X){
  
  X %>%
    as_tibble() %>%
    summarise_all(sum) %>%
    as.matrix()
  
}


p_j <- map(X_j, calculate_p)


m <- 75
rho <- .8


trace_product <- function(A, B) {
  
  a_vec <- as.vector(t(A))
  b_vec <- as.vector(B)
  sum(a_vec * b_vec)
  
}




