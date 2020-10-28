#------------------------------------------------------
# Model-fitting/estimation/testing functions
#------------------------------------------------------

# split a matrix
mat_split <- function(mat, f) {
  lapply(unique(f), function(l) mat[f==l,,drop=FALSE])
}

# trace product 
trace_product <- function(A, B) {
  
  a_vec <- as.vector(t(A))
  b_vec <- as.vector(B)
  sum(a_vec * b_vec)
  
}


robu_handmade <- function(X, y, v, cluster, rho = .8, calc_vcov = NULL) {
  
  k_j <- as.numeric(table(cluster))
  sigma_sq_j <- tapply(v, cluster, mean)
  
  m <- length(k_j)
  
  w_tilde <- 1 / (k_j * sigma_sq_j)
  w_tilde_j <- rep(w_tilde, k_j)
  
  mod_prelim <- lm.wfit(x = X, y = y, w = w_tilde_j)
  
  res <- residuals(mod_prelim) 
  
  # calculate weighted residual sum of squares
  QE <- sum(w_tilde_j * res^2)
  
  # split the design matrix by study
  X_j <- mat_split(X, cluster)
  
  # Create M tilde ----------------------------------------------------------
  M_tilde <- chol2inv(chol(crossprod(X, w_tilde_j * X)))
  p_j <- lapply(X_j, colSums)
  
  # trace products ----------------------------------------------------------
  
  # the first B for the trace product numerator
  w_over_k <- w_tilde / k_j
  B_num_all_1 <- Map(function(w, X) w * crossprod(X), w = w_over_k, X = X_j)
  B_num_1 <- reduce(B_num_all_1, `+`)
  
  # the second B for the trace product numerator
  XJX_j <- lapply(p_j, tcrossprod)
  B_num_all_2 <- Map(function(w, xjx) w * xjx, w = w_over_k, xjx = XJX_j)
  B_num_2 <- reduce(B_num_all_2, `+`)
  
  # the B for the trade product denominator
  B_den_all <- Map(function(w, xjx) w^2 * xjx, w = w_tilde, xjx = XJX_j)
  B_den <- reduce(B_den_all, `+`)
  
  num_minus <- m - (1 - rho) * trace_product(M_tilde, B_num_1) - rho * trace_product(M_tilde, B_num_2)
  den <- sum(k_j * w_tilde) - trace_product(M_tilde, B_den)
  
  tau_sq <- (QE - num_minus) / den
  
  # CE weights
  w_j <- 1 / (k_j * (sigma_sq_j + tau_sq))
  w_ij <- rep(w_j, k_j)
  
  # fit WLS regression
  mod_CE <- lm.wfit(x = X, y = y, w = w_ij)
  
  res <- mod_CE[c("coefficients","residuals","fitted.values","weights")]
  res$tau_sq <- tau_sq
  
  res$X <- X
  res$w_ij <- w_ij
  res$cluster <- cluster
  res$k_j <- k_j
  res$sigma_sq_j <- sigma_sq_j
  res$w_tilde_j <- w_tilde_j
  res$num_minus <- num_minus
  res$den <- den
  res$nobs <- sum(k_j)
  
  class(res) <- "handmade"
  
  if (!is.null(calc_vcov)) {
    res$vcov <- vcovCR(res, cluster = cluster, type = calc_vcov, inverse_var = TRUE)
  }
  
  return(res)
  
}

model.matrix.handmade <- function(object, ...) object$X

bread.handmade <- function(x, ...) {
  x$nobs * chol2inv(chol(crossprod(x$X, x$w_ij * x$X)))
}

update_robu <- function(mod, y, fix_tau_sq = FALSE, calc_CR0 = TRUE) {
  
  if (fix_tau_sq) {
    tau_sq <- mod$tau_sq
  } else {
    mod_prelim <- lm.wfit(x = mod$X, y = y, w = mod$w_tilde_j)
    res <- residuals(mod_prelim) 
    QE <- sum(mod$w_tilde_j * res^2)
    tau_sq <- (QE - mod$num_minus) / mod$den
  }
  
  # CE weights
  w_j <- 1 / (mod$k_j * (mod$sigma_sq_j + tau_sq))
  w_ij <- rep(w_j, mod$k_j)
  
  # fit WLS regression
  mod_CE <- lm.wfit(x = mod$X, y = y, w = w_ij)
  
  res <- mod_CE[c("coefficients","residuals","fitted.values","weights")]
  res$tau_sq <- tau_sq
  
  if (calc_CR0) {
    res$vcov <- calc_CR0(X = mod$X, w_ij = w_ij, 
                         resid = res$residuals, 
                         cluster = mod$cluster)
  }
  
  return(res)
}



# vcovCR ------------------------------------------------------------------

calc_CR0 <- function(X, w_ij, resid, cluster) {
  WX <- w_ij * X
  bread <- chol2inv(chol(crossprod(X, WX)))
  estfun <- apply(resid * WX, 2, function(x) tapply(x, INDEX = cluster, FUN = sum))
  meat <- crossprod(estfun)  
  bread %*% meat %*% bread
}


# run the cwb -------------------------------------------------------------

change_to_mat <- function(res) {
  
  as.matrix(res)
  
}

mult_mat <- function(x, y) {
  
  as.vector(x %*% y)
  
}

mult_mat_const <- function(x, y) {
  
  as.matrix(x * y)
  
}

multiply_rve <- function(X_j, W_j, e_j){
  
  t(X_j) %*% W_j %*% e_j %*% t(e_j) %*% W_j %*% X_j
  
}


calculate_F <- function(beta, vcov, constraints, p = 6, test){
  
  C_mat <- diag(1L, nrow = p)[constraints,,drop = FALSE]    
  
  Q <- as.numeric(t(C_mat %*% beta) %*% chol2inv(chol(C_mat %*% vcov %*% t(C_mat))) %*% C_mat %*% beta)
  q <- nrow(C_mat)
  
  Fstat <- Q/q
  
  test_res <- tibble(test = test, 
                     Fstat = Fstat)
  
  return(test_res)
  
}


cwb <- function(null_model, 
                indices_test, 
                R, 
                full_form, 
                dat) {
  
  
    y <- dat$g
    v <- dat$var_g
    
    X_null <- model.matrix(as.formula(null_model), data = dat)
    X_full <- model.matrix(as.formula(paste("g ~ ", full_form)), data = dat)
    
    cluster <- dat$study
    
    null_mod <- robu_handmade(X = X_null, y = y, v = v, cluster = cluster)
    
    full_mod_org <- robu_handmade(X = X_full, y = y, v = v, cluster = cluster, calc_vcov = "CR0")
    cov_mat_org <- full_mod_org$vcov
    
    # residuals and transformed residuals -------------------------------------
    
    dat$res <- null_mod$residuals
    dat$pred <- null_mod$fitted.values
    split_res <- split(dat$res, dat$study)
    e_tilde_j <- map(split_res, change_to_mat)
    B_j <- attr(vcovCR(null_mod, cluster = cluster, type = "CR2", inverse_var = TRUE), "adjustments")
    dat$t_res <- unlist(pmap(list(B_j, e_tilde_j), mult_mat))
    
    
    # Rademacher weights ------------------------------------------------------
    
    num_studies <- unique(dat$study)
    k_j <- as.numeric(table(dat$study))
    
    bootstraps <- rerun(.n = R, {
      
      wts <- sample(c(-1, 1), size = length(num_studies), replace = TRUE)
      dat$eta <- rep(wts, k_j)
      new_t <- with(dat, pred + res * eta)
      new_t_adj <- with(dat, pred + t_res * eta)
      
      full_mod_cwb <- update_robu(full_mod_org, y = new_t)
      full_mod_cwb_adj <- update_robu(full_mod_org, y = new_t_adj)
      
      cov_mat_cwb <- full_mod_cwb$vcov
      cov_mat_cwb_adj <- full_mod_cwb_adj$vcov
      
      res <- calculate_F(beta = full_mod_cwb$coefficients, vcov = cov_mat_cwb, constraints = indices_test, test = "CWB")
      res_adj <- calculate_F(beta = full_mod_cwb_adj$coefficients, vcov = cov_mat_cwb_adj, constraints = indices_test, test = "CWB Adjusted")
      
      bind_rows(res, res_adj)
      
    }) %>%
      bind_rows()
    
    
    
    org_F <- calculate_F(beta = full_mod_org$coefficients, 
                         vcov = cov_mat_org, 
                         constraints = indices_test, 
                         test = "Naive") %>%
      pull(Fstat)
  
  
  p_boot <- 
    bootstraps %>%
    group_by(test) %>%
    summarize(p_val = mean(Fstat > org_F)) %>%
    ungroup() %>%
    spread(test, p_val)
  
  
  return(p_boot)
  
}


