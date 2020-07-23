library(tidyverse)
library(robumeta)
library(clubSandwich)
library(janitor)
library(broom)
library(estimatr)


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




# Calculate QE ------------------------------------------------------------

# fit model with fixed effects weights
mod_prelim <- lm(delta ~ dv + g2age, data = tsl_dat_small, weights = w_tilde_j)

# residuals
tsl_dat_small <- tsl_dat_small %>%
  mutate(res = residuals(mod_prelim))

# calcualte weighted res sum of squares
QE <- with(tsl_dat_small, sum(w_tilde_j * res^2))

# Create X_j -------------------------------------------------------------

# get the model matrix
full_x <- as.data.frame(model.matrix(mod_prelim))

split(full_x, tsl_dat_small$studyid)

# split the data by studyid
X_j_mat <- full_x %>%
  as_tibble() %>%
  mutate(studyid = tsl_dat_small$studyid) %>%
  group_by(studyid) %>%
  clean_names() %>%
  group_split()


# change to matrix format
change_to_mat <- function(dat){
  
  dat %>%
    select(-studyid) %>%
    as.matrix()
  
}

X_j <- map(X_j_mat, change_to_mat)


# Create W_tilde_j -----------------------------------------------------

# split the w_tilde_j
w_tilde_j <- tsl_dat_small %>%
  select(studyid, w_tilde_j) %>%
  distinct(.) %>%
  group_by(studyid) %>%
  group_split()

w_tilde_j_all <- unique(tsl_dat_small[, c("studyid", "w_tilde_j")])

w_tilde_j <- split(w_tilde_j_all, w_tilde_j_all$studyid)

# make w_tilde_j into vectors

change_to_vec <- function(dat){
  
  dat %>%
    pull(2)
}        



w_tilde_j <- map(w_tilde_j, change_to_vec)


# create identity matrix 
num_es <- tsl_dat_small %>%
  select(studyid, k_j) %>%
  distinct(.) %>%
  select(k_j)

create_identity <- function(k){
  
  diag(k)
  
}

k <- num_es %>% pull(k_j)
I_j <- map(k, create_identity)


# create the big W_j

create_big_W <- function(w, id){
  
  w * id
  
}

W_tilde_j <- pmap(list(w_tilde_j, I_j), create_big_W)



# Create M tilde ----------------------------------------------------------

create_m_tilde <- function(X, W){
  
  t(X) %*% W %*% X
  
}

M_tilde_all <- pmap(list(X_j, W_tilde_j), create_m_tilde)

M_tilde <- M_tilde_all %>% 
  reduce(`+`) %>%
  solve()


# Calculate p_j -----------------------------------------------------------

calculate_p <- function(X){
  
  X %>%
    as_tibble() %>%
    summarise_all(sum) %>%
    as.matrix()
  
}


p_j <- map(X_j, calculate_p)



# trace product -----------------------------------------------------------

weights_tau <- tsl_dat_small %>%
  select(studyid, k_j, w_tilde_j) %>%
  distinct(.) %>%
  mutate(w = w_tilde_j/k_j,
         w_tilde_j_sq = w_tilde_j^2)


wts_num <- weights_tau %>%
  select(studyid, w) %>%
  distinct(.) %>%
  group_by(studyid) %>%
  group_split()

wts_num_j <- map(wts_num, change_to_vec)

wts_den <- weights_tau %>%
  select(studyid, w_tilde_j_sq) %>%
  distinct(.) %>%
  group_by(studyid) %>%
  group_split()

wts_den_j <- map(wts_den, change_to_vec)


calculate_B <- function(w, m){
  
  m_mult <- t(m) %*% m
  w * m_mult
}

# the first B for the trace product numerator
B_num_all_1 <- pmap(list(wts_num_j, X_j), calculate_B)
B_num_1 <- B_num_all_1 %>%
  reduce(`+`)

# the second B for the trace product numerator
B_num_all_2 <- pmap(list(wts_num_j, p_j), calculate_B)
B_num_2 <- B_num_all_2 %>%
  reduce(`+`)

# the B for the trade product denominator
B_den_all <- pmap(list(wts_den_j, p_j), calculate_B)
B_den <- B_den_all %>%
  reduce(`+`)




m <- 119
rho <- .8


trace_product <- function(A, B) {
  
  a_vec <- as.vector(t(A))
  b_vec <- as.vector(B)
  sum(a_vec * b_vec)
  
}

num_tp_1 <- trace_product(M_tilde, B_num_1)
num_tp_2 <- trace_product(M_tilde, B_num_2)
den_tp <- trace_product(M_tilde, B_den)


num <- QE - m + (1 - rho) * num_tp_1 + rho * num_tp_2


den_sum <- weights_tau %>% 
  mutate(den_sum = k_j * w_tilde_j) %>%
  summarize(den_sum = sum(den_sum)) %>%
  pull(den_sum)

den <- den_sum - den_tp

tau_sq <- num/den
tau_sq


# Compare to robu ---------------------------------------------------------

# doesn't match exactly
system.time(robu_comp <- robu(delta ~ g2age + dv, 
                 studynum = studyid, 
                 var.eff.size = v,
                 small = FALSE,
                 data = tsl_dat_small))

dat <- robu_comp$data.full


# estimate  ---------------------------------------------------------------

tsl_dat_small <- tsl_dat_small %>%
  mutate(tau_sq = tau_sq) %>%
  mutate(ce_wts = 1/ (k_j * (sigma_sq_j + tau_sq)))

system.time(mod_full <- lm(delta ~ g2age + dv, data = tsl_dat_small, weights = ce_wts))
mod_null <- lm(delta ~ g2age, data = tsl_dat_small, weights = ce_wts)

mod_res <- tidy(mod_full)

Wald_test(mod_full, constraints = 3:7, vcov = "CR1", cluster = tsl_dat_small$studyid, test = "Naive-F")
Wald_test(mod_full, constraints = 3:7, vcov = "CR2", cluster = tsl_dat_small$studyid, test = "HTZ")

# the HTZ test is a little bit different? 
# HTZ df estimation from lm vs robu is very different 
# because of Var(dst?)
Wald_test(robu_comp, constraints = 3:7, vcov = "CR1", test = "Naive-F")
Wald_test(robu_comp, constraints = 3:7, vcov = "CR2", test = "HTZ")



