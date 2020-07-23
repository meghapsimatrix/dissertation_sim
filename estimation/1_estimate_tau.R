
# change to matrix
change_to_mat <- function(dat){
  
  as.matrix(dat)
  
}

# change to vector
change_to_vec <- function(dat){
  
  as.vector(dat[, 2])
  
}     


# create identity matrix
create_identity <- function(k){
  
  diag(k)
  
}

# create_big_W
create_big_W <- function(w, id){
  
  w * id
  
}

# create M tilde
create_M_tilde <- function(X, W){
  
  t(X) %*% W %*% X
  
}

# calculate p
calculate_p <- function(X){
  
  X %>%
    as_tibble() %>%
    summarise_all(sum) %>%
    as.matrix()
  
}

# calculate B
calculate_B <- function(w, m){
  
  m_mult <- t(m) %*% m
  w * m_mult
}

# trace product 
trace_product <- function(A, B) {
  
  a_vec <- as.vector(t(A))
  b_vec <- as.vector(B)
  sum(a_vec * b_vec)
  
}


# fit model with fixed effects weights

calculate_tau_sq <- function(dat, equation_full, rho = .8){
  
  k_j <- as.numeric(table(dat$study))
  sigma_sq_j <- tapply(dat$v, dat$study, mean)
  w_tilde <- 1 / (k_j * sigma_sq_j)
  dat$w_tilde_j <- w_tilde[dat$study]
  
  m <- length(unique(dat$study))
  
  mod_prelim <- lm(as.formula(equation_full), data = dat, weights = w_tilde_j)
  
  dat$res <- residuals(mod_prelim) 
  
  # calcualte weighted res sum of squares
  QE <- with(dat, sum(w_tilde_j * res^2))
  
  # get the model matrix
  full_x <- as.data.frame(model.matrix(mod_prelim))
  
  # split the data by study
  X_j_mat <- split(full_x, dat$study)
  
  X_j <- map(X_j_mat, change_to_mat)
  
  # Create W_tilde_j -----------------------------------------------------
  
  # split the w_tilde_j
  w_tilde_j_all <- unique(dat[, c("study", "w_tilde_j")])
  w_tilde_j <- split(w_tilde_j_all, w_tilde_j_all$studyid)

  # make w_tilde into vectors
  w_tilde_j <- map(w_tilde_j, change_to_vec)
  
  # create identity matrix 
  num_es <- unique(dat[, c("study", "k_j")])
  k <- as.vector(num_es[, "k_j"])
  I_j <- map(k, create_identity)
  
  
  # create the big W_j
  W_tilde_j <- pmap(list(w_tilde_j, I_j), create_big_W)
  
  # Create M tilde ----------------------------------------------------------
  
  M_tilde_all <- pmap(list(X_j, W_tilde_j), create_M_tilde)
  
  M_tilde <- M_tilde_all %>% 
    reduce(`+`) %>%
    solve()
  
  # Calculate p_j -----------------------------------------------------------
  p_j <- map(X_j, calculate_p)
  
  # trace product -----------------------------------------------------------
  
  weights_tau <- unique(dat[, c("study", "k_j", "w_tilde_j")])
  weights_tau$w <- with(weights_tau, w_tilde_j/k_j)
  weights_tau$w_tilde_j_sq <- with(weights_tau, w_tilde_j^2)
  
  
  wts_num <- unique(weights_tau[, c("study", "w")])
  wts_num <- split(wts_num, weights_tau$studyid)
  
  wts_num_j <- map(wts_num, change_to_vec)
  
  
  wts_den <- unique(weights_tau[, c("study", "w_tilde_j_sq")])
  wts_den <- split(wts_den, weights_tau$studyid)
  
  wts_den_j <- map(wts_den, change_to_vec)
  
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
  

  
  num_tp_1 <- trace_product(M_tilde, B_num_1)
  num_tp_2 <- trace_product(M_tilde, B_num_2)
  den_tp <- trace_product(M_tilde, B_den)
  
  
  num <- QE - m + (1 - rho) * num_tp_1 + rho * num_tp_2
  
  
  den_sum_all <- with(weights_tau, k_j * w_tilde_j)
  den_sum <- sum(den_sum_all)
  
  den <- den_sum - den_tp
  
  tau_sq <- num/den
  
  return(tau_sq)

}
  
  