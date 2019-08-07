pos <- function(X) {
  return(pmax(X, 0))
}

Norm <- function(X, p = 2, invw = 1) {
  # Lp norm
  # invw can be scalar or k-dimensional inverse weight
  # X should be matrix object
  
  n <- nrow(X)
  k <- ncol(X)

  invw <- matrix(rep(invw, each = n), nrow = n, ncol = k)
  Xw = X / invw
 
  return(rowSums(abs(Xw)^p)^(1/p))
}

Vplus <- function(X, kplus_ind) {

  Xplus <- X
  Xplus[,kplus_ind] <- pos(X[, kplus_ind])

  return(Xplus)
}

Vminus <- function(X, kplus_ind) {
  return(Vplus(-X, kplus_ind))
}

maxminQ <- function(mu, sigma, alpha, simlen, lower){
  
  simdata <- MASS::mvrnorm(n = simlen, mu, sigma)
  
  if (lower == T) {
    maxvec <- apply(simdata, 1, max)
    res <- stats::quantile(maxvec, 1-alpha)
  } else {
    minvec <- apply(simdata, 1, min)
    res <- stats::quantile(minvec, alpha)
  }

  return(res)
  
}

maxminQ2 <- function(mu, sigma, alpha, simlen, tau){
  
  simdata <- MASS::mvrnorm(n = simlen, mu, sigma)
  maxvec <- apply(simdata, 1, max)
  res <- stats::quantile(maxvec, 1 - alpha)
  
  return(res - stats::qnorm(1 - tau))
}

gridFun <- function(min_C, max_C, num_grid, J_max, k){
  
  Cgrid <- seq(from = min_C, to = max_C, length.out = num_grid)
  expn <- k / (2 + k)
  C_adj <- seq(from = min_C^expn, to = max_C^expn, length.out = J_max)
  C <- C_adj^(1 / expn)   # Vector of C_j's 
  
  res <- list(Cgrid = Cgrid, C = C)
  return(res)
}

