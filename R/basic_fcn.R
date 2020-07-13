#' Short cut code for max(X,0)
#' 
#' Calculates component-wise max(X,0), where x can be a matrix.
#'
#' @param X A data matrix
#'
#' @return A matrix in the same dimension of X, after pmax(x,0) was applied to each row x of X.
#' @export
#'
#' @examples
pos <- function(X) {
  return(pmax(X, 0))
}

#' Function for (x)_(V+)
#'
#' @param X A data matrix
#' @param mon_ind Monotone variable index, whose length is equal to ncol(X)
#'
#' @return A matrix in the same dimension of X, after (x)_(V+) was applied to each row x of X.
#' @export
#'
#' @examples
Vplus <- function(X, mon_ind) {
  
  Xplus <- X
  Xplus[, mon_ind] <- pos(X[, mon_ind])
  
  return(Xplus)
}

#' Function for (x)_(V-)
#'
#' @param X A numeric matrix
#' @param mon_ind Monotone variable index, whose length is equal to ncol(X)
#'
#' @return A matrix in the same dimension of X, after (x)_(V-) was applied to each row x of X.
#' @export
#'
#' @examples
Vminus <- function(X, mon_ind) {
  
  return(Vplus(-X, mon_ind))
}

#' Function for L_p norm 
#'
#' @param X A numeric matrix (vector is not allowed)
#' @param p Order of the norm
#' @param invw Inverse weights for each component, with the length equal to either 1 or ncol(X)
#'
#' @return A vector of the dimension equal to the nrow(X)
#' @export
#'
#' @examples
Norm <- function(X, p = 2, invw = 1) {
  
  n <- nrow(X)
  d <- ncol(X)
  
  invw <- matrix(rep(invw, each = n), nrow = n, ncol = d)
  Xw <- X / invw
  
  return(rowSums(abs(Xw)^p)^(1/p))
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

