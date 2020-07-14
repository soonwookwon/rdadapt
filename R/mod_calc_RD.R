##' Inverse modulus
##'
##' Calculates the inverse modulus for the regression function at a point
##' problem. More specifically, this calcultes
##' \eqn{\omega^{-1}(b, \Lambda_{\mathcal{V}+}(\gamma_1, C_1),\Lambda_{\mathcal{V}+}(\gamma_2, C_2)) }
##' 
##' @param b point where the inverse modulus is evaluated at
##' @param gam_pair length two vector of gammas (\eqn{(\gamma_1, \gamma_2)'})
##' @param C_pair length two vector of Cs (\eqn{(C_1, C_2)'})
##' @param X n\eqn{\times}k design matrix
##' @param mon_ind index of the monotone variables
##' @param sigma standard deviation of the error term (either length 1 or n)
##' @param swap Indicator for whether we take (C', C) instead of (C, C')
##' 
##' @return A scalar value
##' @export

invmod <- function(b, gam_pair, C_pair, X, mon_ind, sigma = 1, swap = FALSE){
  
  if (!(length(sigma) %in% c(1, nrow(X)))) {
    stop("sigma must have length 1 or n")
  }
  
  K <- b * K_fun(b, gam_pair, C_pair, X, mon_ind, swap)
  res <- sqrt(sum(K^2 / sigma^2))
  
  return(res)
}

##' Calculate omega(0)
##'
##' Calculates the smallest possible b value, or omega(0)
##' 
##' @param gam_pair Not to be used in RD application
##' @param C_pair (C, C')
##' @param X A data matrix
##' @param mon_ind Monotone variable index, whose length is equal to the column length of X 
##' @param swap Indicator for whether we take (C', C) instead of (C, C')

minb_fun <- function(gam_pair, C_pair, X, mon_ind, swap = FALSE){
  
  if (swap) {
    gam_pair <- gam_pair[2:1]
    C_pair <- C_pair[2:1]
  }
  
  minb <- min(C_pair[1] * Norm(Vplus(X, mon_ind))^gam_pair[1] +
                C_pair[2] * Norm(Vminus(X, mon_ind))^gam_pair[2])
  
  return(minb)
}

##' Invsere modulus for RDD
##'
##' Calculates the inverse modulus for the RDD
##' problem. 
##' 
##' @param b point where the inverse modulus is evaluated at
##' @param gam_pair length two vector of gam_pairs (\eqn{(\gam_pair_1, \gam_pair_2)'})
##' @param C_pair length two vector of Cs (\eqn{(C_1, C_2)'})
##' @param Xt \eqn{n_t \times k} design matrix for the treated units
##' @param Xc \eqn{n_c \times k} design matrix for the control units
##' @param mon_ind index of the monotone variables
##' @param sigma_t standard deviation of the error term for the treated units
##' (either length 1 or \eqn{n_t})
##' @param sigma_c standard deviation of the error term for the control units
##' ##' (either length 1 or \eqn{n_c})
##' @param swap indicator of whether to swap the parameter spaces
##' 
##' @return
##' @export

invmod_RD <- function(b, gam_pair, C_pair, Xt, Xc, mon_ind, sigma_t = 1, sigma_c = 1,
                      swap = FALSE){
  
  if (swap) {
    gam_pair <- gam_pair[2:1]
    C_pair <- C_pair[2:1]
  }
  
  ## Derivative of the square of the minization problem 
  deriv_bt <- function(bt) {
    
    bc <- b - bt
    
    om_inv_t_der <- bt * K_fun(bt, gam_pair, C_pair, X, mon_ind) / sigma_t^2
    om_inv_t_der <- sum(om_inv_t_der)
    
    om_inv_c_der <- bc * K_fun(bc, gam_pair, C_pair, X, mon_ind, swap = TRUE) / sigma_c^2
    om_inv_c_der <- sum(om_inv_c_der)
    
    return(om_inv_t_der - om_inv_c_der)
  }
  
  minbt <- minb_fun(gam_pair, C_pair, Xt, mon_ind)
  
  if(b == minbt){
    
    bt <- minbt
  
  }else{
    
    bt_sol <- stats::uniroot(deriv_bt, c(minbt, b), tol = .Machine$double.eps^10)
    bt <- bt_sol$root
  }
  
  delta_t <- invmod(bt, gam_pair, C_pair, Xt, mon_ind, sigma_t)
  delta_c <- invmod(b - bt, gam_pair, C_pair, Xc, mon_ind, sigma_c, swap = TRUE)
  
  res <- list(bt = bt, delta_t = delta_t, bc = b - bt, delta_c = delta_c)
  
  return(res)
  
}

#' Solves for the modulus for the RD parameter
#'
#' @param delta a nonnegative value
#' @param gam_pair not to be used in RD application
#' @param C_pair (C, C')
#' @param Xt \eqn{n_t \times k} design matrix for the treated units
#' @param Xc \eqn{n_c \times k} design matrix for the control units
#' @param mon_ind index of the monotone variables
#' @param sigma_t standard deviation of the error term for the treated units
#' (either length 1 or \eqn{n_t})
#' @param sigma_c standard deviation of the error term for the control units
#' either length 1 or \eqn{n_c})
#' @param swap indicator of whether to swap the parameter spaces
#'
#' @return
#' @export
#'
#' @examples
modsol_RD <- function(delta, gam_pair, C_pair, Xt, Xc, mon_ind, sigma_t, sigma_c,
                      swap = FALSE){
  
  if (swap) {
    gam_pair <- gam_pair[2:1]
    C_pair <- C_pair[2:1]
  }
  
  maxint <- 100 # An arbitrary large number; doesn't affect the result
  
  minbt <- minb_fun(gam_pair, C_pair, Xt, mon_ind)
  minbc <- minb_fun(gam_pair, C_pair, Xc, mon_ind, swap = TRUE)
  minb <- minbt + minbc
  
  eqn_fun <- function(b) {
    
    delta_t <- invmod_RD(b, gam_pair, C_pair, Xt, Xc, mon_ind, sigma_t, sigma_c)$delta_t
    delta_c <- invmod_RD(b, gam_pair, C_pair, Xt, Xc, mon_ind, sigma_t, sigma_c)$delta_c
    
    res <- sqrt(delta_t^2 + delta_c^2) - delta
    return(res)
  }
  
  solve <- stats::uniroot(eqn_fun, c(minb, maxint), extendInt = "upX",
                          tol = .Machine$double.eps^10)
  bsol <- solve$root
  
  res <- invmod_RD(bsol, gam_pair, C_pair, Xt, Xc, mon_ind, sigma_t, sigma_c)
  
  return(res)
}