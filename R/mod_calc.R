##' "Kernel" function for Kwon and Kwon (2019)
##'
##' This function calculates the kernel function that correponds to
##' \eqn{\omega^{-1}(b, \Lambda_{\mathcal{V}+}(\gamma_1, C_1),\Lambda_{\mathcal{V}+}(\gamma_2, C_2)) }
##' 
##' @param b point where the inverse modulus is evaluated at
##' @param gamma length two vector of gammas (\eqn{(\gamma_1, \gamma_2)'})
##' @param C length two vector of Cs (\eqn{(C_1, C_2)'})
##' @param X n\eqn{\times}k design matrix, the kernel is evaluated at each X[i, ]
##' @param mon_ind indice of the monotone variables
##' @param swap indicator of whether to swap the \eqn{\gamma} and \eqn{C}.

K_fun <- function(b, gamma, C, X, mon_ind, swap = FALSE){

  if (swap) {
    gamma <- gamma[2:1]
    C <- C[2:1]
  }
  
  K <- pos(1 - (C[1] / b) * Norm(Vplus(X, mon_ind))^gamma[1] -
             (C[2] / b) * Norm(Vminus(X, mon_ind))^gamma[2])
  
  return(K)
}

##' Calculate \eqn{f_2^* - f_1^*}
##'
##' This function calculates \eqn{f_2^* - f_1^*} that correponds to
##' \eqn{\omega^{-1}(b, \Lambda_{\mathcal{V}+}(\gamma_1, C_1),\Lambda_{\mathcal{V}+}(\gamma_2, C_2)) }
##' 
##' @param b point where the inverse modulus is evaluated at
##' @param gamma length two vector of gammas (\eqn{(\gamma_1, \gamma_2)'})
##' @param C length two vector of Cs (\eqn{(C_1, C_2)'})
##' @param X n\eqn{\times}k design matrix, the kernel is evaluated at each X[i, ]
##' @param mon_ind indice of the monotone variables
##' @param swap indicator of whether to swap the \eqn{\gamma} and \eqn{C}.

f2subf1_fun <- function(b, gamma, C, X, mon_ind, swap = FALSE){
  
  if (swap) {
    gamma <- gamma[2:1]
    C <- C[2:1]
  }
  
  K <- pos(b - (C[1]) * Norm(Vplus(X, mon_ind))^gamma[1] -
             (C[2]) * Norm(Vminus(X, mon_ind))^gamma[2])
  
  return(K)
}

##' Calculate omega(0)
##'
##' Calculates the smallest possible b value, or omega(0)
##' 
##' @param gamma 
##' @param C 
##' @param X 
##' @param mon_ind 
##' @param swap 

minb_fun <- function(gamma, C, X, mon_ind, swap = FALSE){
  
  if (swap) {
    gamma <- gamma[2:1]
    C <- C[2:1]
  }
  
  minb <- min(C[1] * Norm(Vplus(X, mon_ind))^gamma[1] +
                C[2] * Norm(Vminus(X, mon_ind))^gamma[2])
  
  return(minb)
}



##' Invsere modulus
##'
##' Calculates the inverse modulus for the regression function at a point
##' problem. More specifically, this calcultes
##' \eqn{\omega^{-1}(b, \Lambda_{\mathcal{V}+}(\gamma_1, C_1),\Lambda_{\mathcal{V}+}(\gamma_2, C_2)) }
##' 
##' @param b point where the inverse modulus is evaluated at
##' @param gamma length two vector of gammas (\eqn{(\gamma_1, \gamma_2)'})
##' @param C length two vector of Cs (\eqn{(C_1, C_2)'})
##' @param X n\eqn{\times}k design matrix
##' @param mon_ind indice of the monotone variables
##' @param sigma standard deviation of the error term (either length 1 or n)
##' @export

invmod <- function(b, gamma, C, X, mon_ind, sigma = 1, swap = FALSE){
  
  if (!(length(sigma) %in% c(1, nrow(X)))) {
    stop("sigma must have length 1 or n")
  }
  
  # K <- K_fun(b, gamma, C, X, mon_ind, swap = swap)
  # 
  # res <- b * sqrt(sum(K^2 / sigma^2))
  
  K <- f2subf1_fun(b, gamma, C, X, mon_ind, swap = swap)
  
  res <- sqrt(sum(K^2 / sigma^2))
  
  return(res)
}

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



# Calculates ordered modulus of continuity: \omega(\delta)

modsol <- function(delta, gamma, C, X, mon_ind, sigma = 1, swap = FALSE){
  
  maxint <- 100 # fix this
  
  minb <- minb_fun(gamma, C, X, mon_ind)
  
  fun <- function(b) {
    invmod(b, gamma, C, X, mon_ind, sigma) - delta
  }
  
  solve <- stats::uniroot(fun, c(minb, maxint), extendInt = "upX",
                          tol = .Machine$double.eps^10)
  
  return(solve$root)
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

##' Modulus for adaptive RDD
##' 
##' @param gamvec  
##' @param Cvec 
##' @param gam_min
##' @param C_max
##' @param Xt 
##' @param Xc 
##' @param mon_ind 
##' @param sigma_t 
##' @param sigma_c 
##' @param alpha 
##'
##' @export

mod_del_cal <- function(gamvec, Cvec, gam_min = min(gamvec), C_max = max(Cvec), 
                        Xt, Xc, mon_ind, sigma_t, sigma_c, alpha = .05) {

  # For RD design, given C_max = \bar{C}, this function returns
  # 1) \om_t(del_adj,C_J,C_j), \om_c(del_adj,C_J,C_j), \om_t(del_adj,C_j,C_J), \om_c(del_adj,C_j,C_J) 
  # for j=1,...,J, yielding J*4 dim matrix
  # 2) del_jt^L, del_jc^L, del_jt^U, del_jc^U for j=1,...,J, yielding J*4 dim matrix
  # 3) adjusted alpha's for adaptive lower and upper CI 

  J <- length(gamvec)
  
  if (J == 1){   # Minimax case
    
    del_L <- stats::qnorm(1 - alpha)
    del_U <- stats::qnorm(1 - alpha)
    alnewL <- alpha
    alnewU <- alpha
    
  } else {
    
    alnewL <- AdjAlpha_RD(gamvec, Cvec, gam_min, C_max, Xt, Xc, mon_ind, sigma_t, sigma_c,
                          lower = TRUE, alpha = alpha) 
    alnewU <- AdjAlpha_RD(gamvec, Cvec, gam_min, C_max, Xt, Xc, mon_ind, sigma_t, sigma_c,
                          lower = FALSE, alpha = alpha)
    
    del_L <- stats::qnorm(1 - alnewL) 
    del_U <- stats::qnorm(1 - alnewU)   
  }
  
  
  b_mat <- matrix(0, J, 4) # b_tJj, b_cJj, b_tjJ, b_cjJ 
  delta_mat <- matrix(0, J, 4) # Corresponding del_jt^L, del_jc^L, del_jt^U, del_jc^U
  
  for (j in 1:J) {
    
    gampair_j <- c(gamvec[j], gam_min)
    Cpair_j <- c(Cvec[j], C_max)
    
    res_l <- modsol_RD(del_L, gampair_j, Cpair_j, Xt, Xc, mon_ind,
                       sigma_t, sigma_c, swap = TRUE)
    res_u <- modsol_RD(del_U, gampair_j, Cpair_j, Xt, Xc, mon_ind, sigma_t,
                       sigma_c)
    
    b_mat[j, 1:2] <- c(res_l$bt, res_l$bc)
    delta_mat[j, 1:2] <- c(res_l$delta_t, res_l$delta_c)
    b_mat[j, 3:4] <- c(res_u$bt, res_u$bc)
    delta_mat[j, 3:4] <- c(res_u$delta_t, res_u$delta_c)
  }
  
  res <- list(b_mat = b_mat, delta_mat = delta_mat, alnewL = alnewL,
              alnewU = alnewU)

  return(res)
}


# For RD design, given (C_j,C_J) for some fixed j, this function returns
# 1) \om_t(del_adj,C_J,C_j), \om_c(del_adj,C_J,C_j), \om_t(del_adj,C_j,C_J), \om_c(del_adj,C_j,C_J) 
# yielding 1*4 dim matrix
# 2) del_jt^L, del_jc^L, del_jt^U, del_jc^U yielding 1*4 dim matrix

##' mod_del_cal_orc
##'
##' @param gamma
##' @param C
##' @param Xt
##' @param Xc
##' @param mon_ind
##' @param sigma_t
##' @param sigma_c
##' @param alpha
##' @return a list with two elements, with each being a vector of length four
##'   corresponding to 1) \eqn{(\om_t(del,C_J,C_j), \om_c(del,C_J,C_j),
##'   \om_t(del,C_j,C_J), \om_c(del,C_j,C_J))} and 2) \eqn{(del_jt^L, del_jc^L,
##'   del_jt^U, del_jc^U)}, respectively. \eqn{\delta} is fixed at \eqn{z_{1-\alpha}}
##' 
##' @export

mod_del_cal_orc <- function(gamma, C, Xt, Xc, mon_ind, sigma_t, sigma_c,
                            alpha = .05){
  
  del_L <- stats::qnorm(1 - alpha)
  del_U <- stats::qnorm(1 - alpha)
  
  res_l <- modsol_RD(del_L, gamma, C, Xt, Xc, mon_ind, sigma_t, sigma_c,
                     swap = TRUE)
  res_u <- modsol_RD(del_U, gamma, C, Xt, Xc, mon_ind, sigma_t, sigma_c)
  
  b_vec <- c(res_l$bt, res_l$bc, res_u$bt, res_u$bc)
  delta_vec <- c(res_l$delta_t, res_l$delta_c, res_u$delta_t, res_u$delta_c)

  res <- list(b_vec = b_vec, delta_vec = delta_vec)
  
  return(res)  
}

##' Between modulus
##'
##' Calculates the between modulus for the regression function at a point
##' problem
##' 
##' @param delta 
##' @param gamma 
##' @param C 
##' @param X 
##' @param mon_ind 
##' @param sigma
##' @export
bmodsol <- function(delta, gamma, C, X, mon_ind, sigma = 1) {
  
  omega_12 <- modsol(delta, gamma, C, X, mon_ind, sigma)
  omega_21 <- modsol(delta, gamma, C, X, mon_ind, sigma, swap = TRUE)
  
  return(max(omega_12, omega_21))
}




