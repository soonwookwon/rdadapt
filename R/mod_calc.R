##' Invsere modulus
##'
##' Calculates the inverse modulus for the regression function at a point
##' problem. More specifically, this calcultes
##' \eqn{\omega^{-1}(b, \Lambda_{\mathcal{V}+}(\gamma_2, C_2),\Lambda_{\mathcal{V}+}(\gamma_1, C_1)) }
##' 
##' @param b point where the inverse modulus is evaluated at
##' @param gamma length two vector of gammas (\eqn{(\gamma_1, \gamma_2)'})
##' @param C length two vector of Cs (\eqn{(C_1, C_2)'})
##' @param X n\eqn{\times}k design matrix
##' @param mon_ind indice of the monotone variables
##' @param sigma standard deviation of the error term (either length 1 or n)

invmod <- function(b, gamma, C, X, mon_ind, sigma = 1){

  if (!(length(sigma) %in% c(1, nrow(X)))) {
    stop("sigma must have length 1 or n")
  }
  
  C1 <- C[1]
  C2 <- C[2]
  g1 <- gamma[1]
  g2 <- gamma[2]
  
  res <- pos(b - C1 * Norm(Vplus(X,mon_ind))^g1 -
               C2 * Norm(Vminus(X,mon_ind))^g2)^2 
  
  res <- sqrt(sum(res / sigma^2))

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
##' @param ret_list indicator to whether return a more detailed list of results
##' @export

invmod_RD <- function(b, gamma, C, Xt, Xc, mon_ind, sigma_t = 1, sigma_c = 1,
                     ret_list = FALSE){
  # LATER: fix this so that this returns "gf_ip_iota" and others
  # gam = (gamma_1,gamma_2)
  # C = (C_1,C_2)
  # Xt = n_t*k dimensional matrix
  # Xc = n_c*k dimensional matrix
  # sigma_t can be scalar or n_t-dimensional vector for sigma(X_{t,i})
  # sigma_c can be scalar or n_c-dimensional vector for sigma(X_{c,i})
  # When ret_list = 0, the function returns the value of delta
  # When ret_list = 1, the function returns the value of b_t, del_t, b_c, del_c
  
  C1 <- C[1]
  C2 <- C[2]
  g1 <- gamma[1]
  g2 <- gamma[2]
  
  ## Derivative of the square of the minization problem 
  
  deriv_bt <- function(bt) {
    tsum <- pos(bt - C1 * Norm(Vplus(Xt,mon_ind))^g1 -
                  C2 * Norm(Vminus(Xt,mon_ind))^g2) /sigma_t^2
    lhs <- sum(tsum)
    csum <- pos(b - bt - C2 * Norm(Vplus(Xc,mon_ind))^g2 -
                  C1 * Norm(Vminus(Xc,mon_ind))^g1) / sigma_c^2
    rhs <- sum(csum)
    return (lhs - rhs)
  }
  
  if (b == 0) {
    return(0)
  } else {
    btsol <- stats::uniroot(deriv_bt, c(0, b), tol = .Machine$double.eps^10)
    bt <- btsol$root
    
    dt <- invmod(bt, gamma, C, Xt, mon_ind, sigma_t)
    dc <- invmod(b - bt, gamma[2:1], C[2:1], Xc, mon_ind, sigma_c)
    
    if (!ret_list) {
      res <-  sqrt(dt^2 + dc^2)
    } else if (ret_list) {
      res <- list(bt = bt, dt = dt, bc = b - bt, dc = dc)
    }
    
    
    return(res)
  }
  
}



# Calculates ordered modulus of continuity: \omega(\delta)

modsol <- function(delta, gamma, C, X, mon_ind, sigma = 1){
  
  maxint <- 100
  
  fun <- function(b) {
    invmod(b, gamma, C, X, mon_ind, sigma) - delta
  }
  
  solve <- stats::uniroot(fun, c(0, maxint), extendInt = "upX",
                   tol = .Machine$double.eps^10)
  return(solve$root)
}
##'@export
modsol_RD <- function(delta, gamma, C, Xt, Xc, mon_ind, sigma_t, sigma_c,
                      sol_list = FALSE){
  
  maxint <- 100
  
  fun <- function(b) {
    invmod_RD(b, gamma, C, Xt, Xc, mon_ind, sigma_t, sigma_c) - delta
  }
  
  solve <- stats::uniroot(fun, c(0, maxint), extendInt = "upX",
                   tol = .Machine$double.eps^10)
  bsol <- solve$root
  
  if (!sol_list) {
    res <- bsol
  } else if (!sol_list) {
    res <- invmod_RD(bsol, gamma, C, Xt, Xc, mon_ind, sigma_t, sigma_c,
                     ret_list = TRUE)
  }
  
  return(res)
}

# For RD design, given J, this function returns
# 1) \om_t(del_adj,C_J,C_j), \om_c(del_adj,C_J,C_j), \om_t(del_adj,C_j,C_J), \om_c(del_adj,C_j,C_J) 
# for j=1,...,J, yielding J*4 dim matrix
# 2) del_jt^L, del_jc^L, del_jt^U, del_jc^U for j=1,...,J, yielding J*4 dim matrix
# 3) adjusted alpha's for adaptive lower and upper CI 
##'@export
mod_del_cal <- function(gamvec, Cvec, Xt, Xc, mon_ind, sigma_t, sigma_c,
                        alpha = .05) {
  
  jlen <- length(gamvec)
  if (jlen == 1){   # Minimax case
    
    del_L <- stats::qnorm(1-alpha)
    alnewL <- alpha
    alnewU <- alpha
    
  } else {
    
    alnewL <- AdjAlpha_RD(gamvec, Cvec, Xt, Xc, mon_ind, sigma_t, sigma_c,
                          lower = TRUE, alpha = alpha)   # adjusted alphas
    
    alnewU <- AdjAlpha_RD(gamvec, Cvec, Xt, Xc, mon_ind, sigma_t, sigma_c,
                          lower = FALSE, alpha = alpha)
    del_L <- stats::qnorm(1-alnewL)   # adjusted deltas (corresponds to \delta^{adpt})
    del_U <- stats::qnorm(1-alnewU)   
  }
  
  
  b_mat <- matrix(0, jlen, 4) # b_tJj, b_cJj, b_tjJ, b_cjJ 
  d_mat <- matrix(0, jlen, 4) # Corresponding del_jt^L, del_jc^L, del_jt^U, del_jc^U
  
  for (j in 1:jlen) { 
    gampair_j <- c(gamvec[j], gamvec[jlen])
    Cpair_j <- c(Cvec[j], Cvec[jlen])
    
    res_l <- modsol_RD(del_L, gampair_j[2:1], Cpair_j[2:1], Xt, Xc, mon_ind,
                       sigma_t, sigma_c, sol_list = TRUE)
    
    if(jlen == 1){   # minimax case
      
      res_u <- res_l
      
    }else{
      
      res_u <- modsol_RD(del_U, gampair_j, Cpair_j, Xt, Xc, mon_ind, sigma_t,
                         sigma_c, sol_list = 1)
      
    }
    
    b_mat[j, 1:2] <- c(res_l$bt, res_l$bc)
    d_mat[j, 1:2] <- c(res_l$dt, res_l$dc)
    b_mat[j, 3:4] <- c(res_u$bt, res_u$bc)
    d_mat[j, 3:4] <- c(res_u$dt, res_u$dc)
  }
  
  res <- list(b_mat = b_mat, d_mat = d_mat, alnewL = alnewL, alnewU = alnewU)
  return(res)
}


# For RD design, given (C_j,C_J) for some fixed j, this function returns
# 1) \om_t(del_adj,C_J,C_j), \om_c(del_adj,C_J,C_j), \om_t(del_adj,C_j,C_J), \om_c(del_adj,C_j,C_J) 
# yielding 1*4 dim matrix
# 2) del_jt^L, del_jc^L, del_jt^U, del_jc^U yielding 1*4 dim matrix

mod_del_cal_orc <- function(gamma, C, maxgam, maxC, Xt, Xc, mon_ind,
                            sigma_t, sigma_c, alpha = .05){
  
  del_L <- stats::qnorm(1 - alpha)
  del_U <- stats::qnorm(1 - alpha)
  
  b_mat <- matrix(0, 1, 4) # b_tJj, b_cJj, b_tjJ, b_cjJ 
  d_mat <- matrix(0, 1, 4) # Corresponding del_jt^L, del_jc^L, del_jt^U, del_jc^U
  
  gampair_j <- c(gamma, maxgam)
  Cpair_j <- c(C, maxC)
  
  res_l <- modsol_RD(del_L, gampair_j[2:1], Cpair_j[2:1], Xt, Xc, mon_ind,
                     sigma_t, sigma_c, sol_list = TRUE)
  res_u <- modsol_RD(del_U, gampair_j, Cpair_j, Xt, Xc, mon_ind,
                     sigma_t, sigma_c, sol_list = TRUE)
  
  b_mat[1, 1:2] <- c(res_l$bt, res_l$bc)
  d_mat[1, 1:2] <- c(res_l$dt, res_l$dc)
  b_mat[1, 3:4] <- c(res_u$bt, res_u$bc)
  d_mat[1, 3:4] <- c(res_u$dt, res_u$dc)
  
  res <- list(b_mat = b_mat, d_mat = d_mat)
  return(res)  
}


# Calculates between modulus of continuity: \omega_+(\delta)

bmodsol <- function(delta, gamma, C, X, mon_ind, sigma = 1) {
  
  res1 <- modsol(delta, gamma, C, X, mon_ind, sigma)
  res2 <- modsol(delta, gamma[c(2,1)], C[c(2,1)], X, mon_ind, sigma)
  
  return(max(res1, res2))
}




