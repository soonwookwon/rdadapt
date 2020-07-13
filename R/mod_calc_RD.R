##' Invsere modulus for RDD
##'
##' Calculates the inverse modulus for the RDD
##' problem. 
##' 
##' @param b point where the inverse modulus is evaluated at
##' @param gamma length two vector of gammas (\eqn{(\gamma_1, \gamma_2)'})
##' @param C length two vector of Cs (\eqn{(C_1, C_2)'})
##' @param Xt \eqn{n_t \times k} design matrix for the treated units
##' @param Xc \eqn{n_c \times k} design matrix for the control units
##' @param mon_ind indice of the monotone variables
##' @param sigma_t standard deviation of the error term for the treated units
##' (either length 1 or \eqn{n_t})
##' @param sigma_c standard deviation of the error term for the control units
##' @param swap indicator of whether to swap the parameter spaces
##' 
##' @export

invmod_RD <- function(b, gamma, C, Xt, Xc, mon_ind, sigma_t = 1, sigma_c = 1,
                      swap = FALSE){
  
  if (swap) {
    gamma <- gamma[2:1]
    C <- C[2:1]
  }
  
  ## Derivative of the square of the minization problem 
  deriv_bt <- function(bt) {
    
    om_inv_t_der <- f2subf1_fun(bt, gamma, C, Xt, mon_ind) / sigma_t^2
    om_inv_t_der <- sum(om_inv_t_der)
    
    om_inv_c_der <- f2subf1_fun(b - bt, gamma, C, Xc, mon_ind, swap = TRUE) /
      sigma_c^2
    om_inv_c_der <- sum(om_inv_c_der)
    
    return(om_inv_t_der - om_inv_c_der)
  }
  
  minbt <- minb_fun(gamma, C, Xt, mon_ind)
  
  if(b == minbt){
    bt <- minbt
  }else{
    bt_sol <- stats::uniroot(deriv_bt, c(minbt, b), tol = .Machine$double.eps^10)
    
    bt <- bt_sol$root
  }
  
  delta_t <- invmod(bt, gamma, C, Xt, mon_ind, sigma_t)
  delta_c <- invmod(b - bt, gamma, C, Xc, mon_ind, sigma_c, swap = TRUE)
  
  res <- list(delta = sqrt(delta_t^2 + delta_c^2), bt = bt, delta_t = delta_t,
              bc = b - bt, delta_c = delta_c)
  
  return(res)
  
}


##' Calculate modulus for RDD
##'
##' @param delta 
##' @param gamma 
##' @param C 
##' @param Xt 
##' @param Xc 
##' @param mon_ind 
##' @param sigma_t 
##' @param sigma_c 
##' @param sol_list
##' @export

modsol_RD <- function(delta, gamma, C, Xt, Xc, mon_ind, sigma_t, sigma_c,
                      swap = FALSE){
  
  if (swap) {
    gamma <- gamma[2:1]
    C <- C[2:1]
  }
  
  maxint <- 100
  
  minbt <- minb_fun(gamma, C, Xt, mon_ind)
  minbc <- minb_fun(gamma, C, Xc, mon_ind, swap = TRUE)
  minb <- minbt + minbc
  
  fun <- function(b) {
    invmod_RD(b, gamma, C, Xt, Xc, mon_ind, sigma_t, sigma_c)$delta - delta
  }
  
  solve <- stats::uniroot(fun, c(minb, maxint), extendInt = "upX",
                          tol = .Machine$double.eps^10)
  
  bsol <- solve$root
  
  res <- invmod_RD(bsol, gamma, C, Xt, Xc, mon_ind, sigma_t, sigma_c)
  
  return(res)
}