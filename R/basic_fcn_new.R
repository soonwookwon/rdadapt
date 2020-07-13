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

#' Function for K(x/b; C, C')
#'
#' @param b "Bandwidth" parameter
#' @param gam_pair Not to be used in RDD paper
#' @param C_pair (C, C')
#' @param X A data matrix
#' @param mon_ind Monotone variable index, whose length is equal to the column length of X 
#' @param swap Indicator for whether we take (C', C) instead of (C, C')
#'
#' @return
#' @export
#'
#' @examples
K_fun <- function(b, gam_pair = c(1, 1), C_pair, X, mon_ind, swap = FALSE){
  
  if (swap) {
    gam_pair <- gam_pair[2:1]
    C_pair <- C_pair[2:1]
  }
  
  K <- pos(1 - (C_pair[1] / b) * Norm(Vplus(X, mon_ind))^gam_pair[1] -
             (C_pair[2] / b) * Norm(Vminus(X, mon_ind))^gam_pair[2])
  
  return(K)
}
