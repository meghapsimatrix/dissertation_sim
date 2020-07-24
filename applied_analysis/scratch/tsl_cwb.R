library(tidyverse)
library(robumeta)
library(clubSandwich)


# Load data and clean -----------------------------------------------------

load("data/tsl_dat.RData")

# complete cases 
tsl_dat_small <- tsl_dat %>%
  select(studyid, esid, delta, v, dv, g2age) %>%
  drop_na()

# Full and null model -----------------------------------------------------

cw_full <- robu(delta ~ dv + g2age, 
                studynum = studyid, 
                var.eff.size = v,
                small = TRUE,
                data = tsl_dat_small)

Wald_test(cw_full, constraints = 2:6, vcov = "CR2", test = "Naive-F")
Wald_test(cw_full, constraints = 2:6, vcov = "CR2", test = "HTZ")

cw_null <- robu(delta ~ g2age,
                studynum = studyid, 
                var.eff.size = v,
                data = tsl_dat_small,
                small = TRUE)



# Brute force -------------------------------------------------------------

# residuals from null 
tsl_dat_small$res_null <- clubSandwich:::residuals_CS.robu(cw_null)




# creating the -1 and 1 randomly generated weights for the residuals  per study
num_studies <- unique(tsl_dat_small$studyid)
eta <- sample(c(-1, 1), size = length(num_studies), replace = TRUE)

B_j <- attr(vcovCR(cw_null, type = "CR2"), "adjustments")

eta_dat <- tibble(studyid = num_studies, eta = eta)

tsl_dat_small <- left_join(tsl_dat_small, eta_dat, by = "studyid")

tsl_dat_small <- tsl_dat_small %>%
  mutate(eta = eta, 
         new_y = delta + res_null * eta)  # I need to multiply by B_j somehow

# re-estimating

cw_reestimate <- robu(new_y ~ dv + g2age, 
                      studynum = studyid, 
                      var.eff.size = v,   # do I need to recalculate v??
                      small = TRUE,
                      data = tsl_dat_small)

Wald_test(cw_reestimate, constraints = 2:6, vcov = "CR2", test = "Naive-F")
Wald_test(cw_reestimate, constraints = 2:6, vcov = "CR2", test = "HTZ")

# need to do boostrapping somewhere here



# Computationally efficient -----------------------------------------------

# do the regression stuff with the design matrix 
# same formula 
# one variable at a time 
# use the same weight matrix 

# x matrix
full_x <- clubSandwich:::model_matrix.robu(cw_full)

head(full_x)
str(full_x)



# Create the Uj and Xj matrix  --------------------------------------------


Uj <- full_x[, c(1, 7)] %>%
  as_tibble() %>%
  mutate(studyid = tsl_dat_small$studyid) %>%
  group_by(studyid)

Xj <- full_x[, 2:6] %>%
  as_tibble() %>%
  mutate(studyid = tsl_dat_small$studyid) %>%
  group_by(studyid)

U_j_mat <- group_split(Uj)
X_j_mat <- group_split(Xj)


change_to_mat <- function(dat){
  
  dat %>%
    select(-studyid) %>%
    as.matrix()
  
}


U_j <- map(U_j_mat, change_to_mat)
X_j <- map(X_j_mat, change_to_mat)

head(U_j)
head(X_j)


# Weight matrix -----------------------------------------------------------

# block diagonal 
W_j <- clubSandwich:::weightMatrix.robu(cw_full, cluster = cw_full$study_orig_id)

head(W_j)

is.matrix(W_j[[119]])

# The outcome -------------------------------------------------------------
Tj <- tsl_dat_small %>%
  select(delta, studyid) %>%
  group_by(studyid)

T_j_mat <- Tj %>%
  group_split()

T_j <- map(T_j_mat, change_to_mat)

# Adjustment matrices -----------------------------------------------------
# adjustment matrices
cr2_full <- vcovCR(cw_full, type = "CR2")

# block diagonal - A_j  - from full model 
A_j <- attr(cr2_full, "adjustments")


# Calculation by hand ------------------------------------------------------

# Calculate M_u

calculate_m <- function(var, wt) {
  
  t(var) %*% wt %*% var
  
}

M_u_all <- pmap(list(U_j, W_j), calculate_m)
M_u <- solve(M_u_all %>% reduce(`+`)) 

# calculate a tilde

calculate_coefs <- function(var, weights, outcome) {
  
  t(var) %*% weights %*% outcome
  
}

a_all <- pmap(list(U_j, W_j, T_j), calculate_coefs)
a_tilde <- a_all %>% 
  reduce(`+`)

a_tilde <- M_u %*% a_tilde  # is this right?

# calculate e tilde 

calculate_e <- function(outcome, var, coefs) {
  
  outcome - var %*% coefs
  
}

e_all <- pmap(list(T_j, U_j, a_tilde %>% list()), calculate_e)

e_tilde <- e_all



# X” = (I – H_U) X = X – H_U X

# In other words, calculate 
#H_U X = U M_U U’ W X (the fitted values) and then subtract that from X to get the residuals. 


calculate_hx <- function(var, M, weights, var_2){
  
  var %*% M %*% t(var) %*% weights %*% var_2
  
}


hx_all <- pmap(list(U_j, M_u %>% list(), W_j, X_j), calculate_hx) # ?????



get_res <- function(var, var_2){
  
  var - var_2
  
}

x_dot_dot_j <- pmap(list(X_j, hx_all), get_res) # ?????


# beta hat

M_x_dot_dot_all <- pmap(list(x_dot_dot_j, W_j), calculate_m)
M_x_dot_dot <- solve(M_x_dot_dot_all %>% reduce(`+`)) 


b_all <- pmap(list(x_dot_dot_j, W_j, e_tilde), calculate_coefs)
b_hat <- b_all %>% 
  reduce(`+`)

b_hat <- M_x_dot_dot %*% b_hat


# ehat
e_hat <- pmap(list(e_tilde, x_dot_dot_j, b_hat %>% list()), calculate_e)



# Variance robust calculate by hand 


calculate_vr_meat <- function(var, weights, adj, res){
  
  t(var) %*% weights %*% adj %*% res %*% t(res) %*% adj %*% weights %*% var
  
}

vr_meat_all <- pmap(list(x_dot_dot_j, W_j, A_j, e_hat), calculate_vr_meat)

vr_meat <- vr_meat_all %>% 
  reduce(`+`)

Vr <- M_x_dot_dot %*% vr_meat %*% M_x_dot_dot



# Calculate Q

Q <- t(b_hat) %*% solve(Vr) %*% b_hat


# T star

calculate_t_star <- function(var_1, coef_1, r_wts, adj, res){
  
  var_1 %*% coef_1 + r_wts * adj %*% res
  
}

eta_j <- tsl_dat_small %>%
  select(eta, studyid) %>%
  group_by(studyid)

eta_j_mat <- eta_j %>%
  group_split()

eta_j <- map(eta_j_mat, change_to_mat)


T_star_j <- pmap(list(U_j, a_tilde %>% list(), eta_j, B_j, e_tilde), calculate_t_star)


# f_j


calculate_f <- function(var, res){
  
  var %*% res
  
}

f_j <- pmap(list(B_j, e_tilde), calculate_f)


calculate_e_j <- function(r_wts, f){
  
  r_wts * f
}


e_j_r <- pmap(list(eta_j, f_j), calculate_e_j)


calculate_E_t <- function(var, wts){

 t(var) %*% wts

}


calculate_G_t <- function(var, wts){
  
  var %*% wts
  
}




E_t_j <- pmap(list(x_dot_dot_j, W_j), calculate_E_t)

G_t_j <- pmap(list(E_t_j, A_j), calculate_G_t)


calculate_b_hat_r <- function(var, res){
  
  var %*% res
  
}


# in 2.43 should be tranpose of E because non conformable
b_hat_r_part <- pmap(list(E_t_j, e_j_r), calculate_b_hat_r)

b_hat_r_part <- b_hat_r_part %>% 
  reduce(`+`)

b_hat_r <- M_x_dot_dot %*% b_hat_r_part


e_hat_j_r <- pmap(list(e_j_r, x_dot_dot_j, b_hat_r %>% list()), calculate_e)


calculate_V_bootstap <- function(var, res){
  
  var %*% res %*% t(res) %*% t(var)
  
}


V_bootstrap_meat_all <- pmap(list(G_t_j, e_hat_j_r), calculate_V_bootstap)

V_bootstrap_meat <- V_bootstrap_meat_all %>% 
  reduce(`+`)

V_r <- M_x_dot_dot %*% V_bootstrap_meat %*% M_x_dot_dot



Q_r <- t(b_hat_r) %*% solve(V_r) %*% b_hat_r


# do this for bootstrapping

