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
                small = TRUE,
                data = tsl_dat_small)


# Brute force -------------------------------------------------------------

# residuals from null 
tsl_dat_small$res_null <- clubSandwich:::residuals_CS.robu(cw_null)


# creating the -1 and 1 randomly generated weights for the residuals  per study
num_studies <- unique(tsl_dat_small$studyid)
p <- rbinom(n = length(num_studies), size = 1, prob = .5)

eta_dat <- tibble(studyid = num_studies, p = p) %>%
  mutate(eta = ifelse(p == 1, 1, -1))

tsl_dat_small <- left_join(tsl_dat_small, eta_dat, by = "studyid")

tsl_dat_small <- tsl_dat_small %>%
  mutate(new_y = delta + res_null * eta)  # is this correct??

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

Xj <- full_x[, 1:6] %>%
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
M_u <- solve(M_u_all %>% reduce(`+`)) %>% 
  list()

# calculate a tilde

calculate_a_tilde <- function(M, var, weights, outcome) {
  
  M %*% t(var) %*% weights %*% outcome
  
}

a_all <- pmap(list(M_u, U_j, W_j, T_j), calculate_a_tilde)
a_tilde <- a_all %>% 
  reduce(`+`) %>% 
  list()

# calculate e tilde 

calculate_e_tilde <- function(outcome, var, coefs) {
  
  outcome - var %*% coefs
  
}

calculate_h <- function(var, M, weights){
  
  var %*% M %*% t(var) %*% weights
  
}

# how do I do the identity matrix thing here (I - H_u) T  ???

e_all <- pmap(list(T_j, U_j, a_tilde), calculate_e_tilde)
e_tilde <- e_all %>% 
  unlist()


h_all <- pmap(list(U_j, M_u, W_j), calculate_h) # ?????
H_u <- h_all %>% 
  reduce(`+`) %>% 
  list()

# x dot dot - do I do this per study and then combine somehow??


