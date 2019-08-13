##' Minimax CI length for regression function at a point
##'
##' This function calculates the minimax CI length for the regression function
##' at a point problem
##' 
##' @param b 
##' @param gamma 
##' @param C 
##' @param X 
##' @param mon_ind 
##' @param sigma 
##' @param alpha
##' @export

CI_length <- function(b, gamma, C, X, mon_ind, sigma = 1, alpha = .05) {

  om_inv <- invmod(b, rep(gamma, 2), rep(C, 2), X, mon_ind, sigma)

  if (om_inv == 0) {
    return(Inf)
  }

  K <- K_fun(b, gamma, C, X, mon_ind)
  gf_ip_iota <- sum(b * K / sigma^2)
  
  sd <- om_inv / gf_ip_iota
  bias <- .5 * ( b - ( om_inv^2 / gf_ip_iota ))
  
  cva <- ifelse((bias / sd) > 3,
                (bias / sd) + stats::qnorm(1 - alpha),
                sqrt(stats::qchisq(1 - alpha, df = 1, ncp = (bias / sd)^2)))
  
  return(2 * cva * sd)
}

##' Length of minimax CI for the RDD
##'
##' Calculates the length of the minimax CI for the RDD
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
##' @param alpha significance level
##' @export

CI_length_RD <- function(b, gamma, C, Xt, Xc, mon_ind, sigma_t = 1, sigma_c = 1,
                         alpha = .05) {
  
  om_inv <- invmod_RD(b, rep(gamma,2), rep(C,2), Xt, Xc, mon_ind = mon_ind,
                      sigma_t = sigma_t, sigma_c = sigma_c)

  bc <- om_inv$bc
  bt <- om_inv$bt
  om_inv_t <- om_inv$delta_t 
  om_inv_c <- om_inv$delta_c
  om_inv <- om_inv$delta
    
  if (om_inv == 0) return(Inf)

  Kt <- K_fun(bt, gamma, C, Xt, mon_ind)
  gf_ip_iota_t <- sum(bt * Kt / sigma_t^2)

  Kc <- K_fun(bc, gamma, C, Xc, mon_ind)
  gf_ip_iota_c <- sum(bc * Kc / sigma_c^2)
  
  sd <- sqrt((om_inv_t/gf_ip_iota_t)^2 + (om_inv_c/gf_ip_iota_c)^2)
  
  bias <- .5 * (b - (om_inv^2 /  (gf_ip_iota_t + gf_ip_iota_c)))
  cva <- ifelse(abs(bias / sd) > 3,
                abs(bias / sd) + stats::qnorm(1 - alpha),
                sqrt(stats::qchisq(1 - alpha, df = 1, ncp = (bias / sd)^2)))
  
  return(2 * cva * sd)
}

##' Estimator used in constructing the minimax (or adaptive) CIs 
##'
##' Calculates the estimator used in constructing then minimax (or adaptive) CIs
##' for the regression function at a point problem. Again, the
##' estimator corresponds to \eqn{\omega^{-1}(b,
##' \Lambda_{\mathcal{V}+}(\gamma_1, C_1),\Lambda_{\mathcal{V}+}(\gamma_2, C_2))
##' }
##' 
##' 
##' @param b 
##' @param gamma 
##' @param C 
##' @param X 
##' @param mon_ind 
##' @param sigma 
##' @param Y 
##' @param swap 
##' @export
Lhat_fun <- function(b, gamma, C, X, mon_ind, sigma, Y, swap = FALSE) {

  if (swap) {
    gamma <- gamma[2:1]
    C <- C[2:1]
  } 
  
  f1 <- C[1] * Norm(Vplus(X, mon_ind))^gamma[1]
  f2 <- pmax(b - C[2] * Norm(Vminus(X, mon_ind))^gamma[2], f1)
    
  return(.5 * b + sum((f2 - f1) * (Y - .5 * (f1 + f2)) / sigma^2) /
           sum((f2 - f1) / sigma^2))
}


##' Estimator used in constructing the minimax CI for RDD
##'
##' Calculates the estimator that is used to construct the adaptive (or minimax)
##' CI for RDD. Again, the estimator correponds to the inverse modulus
##' \eqn{\omega^{-1}(b, \Lambda_{\mathcal{V}+}(\gamma_1,
##' C_1),\Lambda_{\mathcal{V}+}(\gamma_2, C_2)). }
##' 
##' @param b point where the inverse modulus is evaluated at
##' @param gamma length two vector of gammas (\eqn{(\gamma_1, \gamma_2)'})
##' @param C length two vector of Cs (\eqn{(C_1, C_2)'})
##' @param Xt \eqn{n_t \times k} design matrix for the treated units
##' @param Xc \eqn{n_c \times k} design matrix for the control units
##' @param Yt length \eqn{n_t} of outcome variables for the treated units
##' @param Yc length \eqn{n_c} of outcome variables for the control units
##' @param mon_ind indice of the monotone variables
##' @param sigma_t standard deviation of the error term for the treated units
##' (either length 1 or \eqn{n_t})
##' @param sigma_c standard deviation of the error term for the control units
##' @param swap indiactor of whether to swap the parameter spaces
##' @export
Lhat_RD_fun <- function(b = NULL, bt = NULL, bc = NULL, gamma, C, Xt, Xc, Yt, Yc,
                        mon_ind, sigma_t, sigma_c, swap = FALSE) {

  if (swap) {
    gamma <- gamma[2:1]
    C <- C[2:1]
  }

  if (is.null(bt) | is.null(bc)) {

    if (is.null(b)) {
      stop("You have to provide either 1) b or 2) bt and bc")
    }
    
    om_inv <- invmod_RD(b, gamma, C, Xt, Xc, mon_ind = mon_ind,
                        sigma_t = sigma_t, sigma_c = sigma_c)
    bt <- om_inv$bt
    bc <- om_inv$bc
  }
  
 
  Lhat_t <- Lhat_fun(bt, gamma, C, Xt, Yt, mon_ind = mon_ind, sigma = sigma_t)
  Lhat_c <- Lhat_fun(bc, gamma, C, Xc, Yc, mon_ind = mon_ind, sigma = sigma_c,
                     swap = TRUE)

  return(Lhat_t - Lhat_c)
}

##' (One-sided) Adaptive CI for RDD
##'
##' Calculates the (one-sided) adaptive CI for RDD 
##' 
##' @param bpairmat 
##' @param dmat 
##' @param gamma 
##' @param C 
##' @param gam_min
##' @param C_max
##' @param Xt 
##' @param Xc 
##' @param mon_ind 
##' @param sigma_t 
##' @param sigma_c 
##' @param Yt 
##' @param Yc 
##' @param lower 
##' @param alpha 
##' @export
CI_one_sd_RD <- function(bpairmat, dmat, gamma, C, gam_min = min(gamma), C_max = max(C), 
                         Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, lower, alpha = .05){
  
  J <- length(gamma)
  
  bpairmat_l <- bpairmat[, 1:2, drop=F]
  bpairmat_u <- bpairmat[, 3:4, drop=F]
  dmat_l <- dmat[, 1:2, drop=F]
  dmat_u <- dmat[, 3:4, drop=F]
  
  chatvec <- numeric(J)
  
  
  for (j in 1:J) {

    gampair_j <- c(gamma[j], gam_min)
    Cpair_j <- c(C[j], C_max)
    
    if (lower) {

      swap <- TRUE
      
      bt <- bpairmat_l[j, 1]
      bc <- bpairmat_l[j, 2]
      delta_t <- dmat_l[j, 1]
      delta_c <- dmat_l[j, 2]
      
    } else {

      swap <- FALSE
      
      bt <- bpairmat_u[j, 1]
      bc <- bpairmat_u[j, 2]
      delta_t <- dmat_u[j, 1]
      delta_c <- dmat_u[j, 2]
      
    } 

    K_jt <- K_fun(bt, gampair_j, Cpair_j, Xt, mon_ind = mon_ind, swap = TRUE)
    K_jc <- K_fun(bc, gampair_j, Cpair_j, Xc, mon_ind = mon_ind)
    
    omega_der_t <- delta_t / sum(bt * K_jt / sigma_t^2) 
    omega_der_c <- delta_c / sum(bc * K_jc / sigma_c^2)
    omega_der <- sqrt(omega_der_t^2 + omega_der_c^2)
    
    Lhat <- Lhat_RD_fun(bt = bt, bc = bc, gamma = gampair_j, C = Cpair_j,
                        Xt = Xt, Xc = Xc, Yt = Yt, Yc = Yc, mon_ind = mon_ind,
                        sigma_t = sigma_t, sigma_c = sigma_c, swap = swap)
    
    chat <- ifelse(lower == TRUE, 
                   Lhat - 0.5 * ((bt + bc) + stats::qnorm(1 - alpha) * omega_der), 
                   Lhat + 0.5 * ((bt + bc) + stats::qnorm(1 - alpha) * omega_der))
    
    
    chatvec[j] <- chat  
  }
  
  res <- ifelse(lower == T, max(chatvec), min(chatvec))
 
  return(res) 
}

##'  Adjusting alpha for one-sided CI for regression function at a point problem
##'
##' @param gamma 
##' @param C 
##' @param X 
##' @param sigma 
##' @param mon_ind 
##' @param lower 
##' @param simlen 
##' @param alpha 
AdjAlpha <- function(gamma, C, X, sigma, mon_ind, lower, simlen = 1e04,
                     alpha = .05){
  
  J <- length(gamma)
  n <- length(X[, 1])
  
  f <- function(a) {
    
    del_adpt <- stats::qnorm(1 - a)
    b_mat_bon <- matrix(0, J, 2)
    
    for(j in 1:J){
      gampair_j <- c(gamma[j], gamma[J])
      Cpair_j <- c(C[j], C[J])
      
      b_mat_bon[j, 1] <- modsol(del_adpt, gampair_j[2:1], Cpair_j[2:1], X,
                                mon_ind, sigma = sigma)
      b_mat_bon[j, 2] <- modsol(del_adpt, gampair_j, Cpair_j, X, mon_ind,
                                sigma = sigma)
    }
    
    sumpartmat <- matrix(0, n, J)
    omprvec <- numeric(J)
    
    delta <- del_adpt
    
    for (j in 1:J) {
      
      if (lower == TRUE) {
        b <- b_mat_bon[j, 1]
        C1 <- C[J]
        C2 <- C[j]
        g1 <- gamma[J]
        g2 <- gamma[j]
      } else {
        b <- b_mat_bon[j, 2]
        C1 <- C[j]
        C2 <- C[J]
        g1 <- gamma[j]
        g2 <- gamma[J]
      } 

      f1_minus_f2 <- pos(b - C1 * Norm(Vplus(X, mon_ind))^g1 - 
                           C2 * Norm(Vminus(X, mon_ind))^g2) / sigma
      omega_der_denom <- sum(f1_minus_f2 / sigma)
      
      omega_der <- delta / omega_der_denom
      sumpartmat[, j] <- f1_minus_f2
      omprvec[j] <- omega_der
      
    }
      
    muvec <- ifelse(lower==T,  stats::qnorm(a),  stats::qnorm(1-a)) * omprvec
    sigmamat <- matrix(0, J, J)
    
    for (j1 in 1:J) {
      for (j2 in j1:J) {
        if (j1 == j2) {
          sigmamat[j1, j1] <- omprvec[j1]^2
        } else {
          sigmamat[j1, j2] <- (omprvec[j1] * omprvec[j2] / delta^2) *
            sum(sumpartmat[, j1] * sumpartmat[, j2])
          sigmamat[j2, j1] <- sigmamat[j1, j2]
        }
      }
    }
      
    q <- maxminQ(muvec, sigmamat, alpha/2, simlen, lower)
    return(q)
  }
    
  r <- stats::uniroot(f, c(alpha / (2 * J), alpha / 2))
  res <- r$root
  
  return(res)
}



##' Adjusting alpha for two-sided CI for regression function at a point problem
##'
##' @param gamma 
##' @param C 
##' @param X 
##' @param sigma 
##' @param mon_ind 
##' @param simlen 
##' @param alpha 
AdjAlpha2 <- function(gamma, C, X, sigma, mon_ind, simlen = 1e04,
                      alpha = .05) {
  
  J <- length(gamma)
  n <- nrow(X)
  
  f <- function(a){
    
    del_adpt <- stats::qnorm(1 - a)
    b_mat_bon <- matrix(0, J, 2)
    
    for(j in 1:J){
      gampair_j <- c(gamma[j], gamma[J])
      Cpair_j <- c(C[j], C[J])
      
      b_mat_bon[j, 1] <- modsol(del_adpt, gampair_j[2:1], Cpair_j[2:1], X,
                                sigma = sigma, mon_ind)
      b_mat_bon[j, 2] <- modsol(del_adpt, gampair_j, Cpair_j, X,
                                sigma = sigma, mon_ind)
    }
    
    sumpartmatL <- matrix(0, n, J)
    omprvecL <- numeric(J)
    sumpartmatU <- matrix(0, n, J)
    omprvecU <- numeric(J)
    
    delta <- del_adpt
    
    for(j in 1:J){
      
      b <- b_mat_bon[j, 1]
      C1 <- C[J]
      C2 <- C[j]
      g1 <- gamma[J]
      g2 <- gamma[j]

      f1_minus_f2 <- pos(b - C1 * Norm(Vplus(X, mon_ind))^g1 - 
                           C2 * Norm(Vminus(X, mon_ind))^g2) / sigma      
      omega_der_denom <- sum(f1_minus_f2 / sigma)
      omega_der <- delta / omega_der_denom
      
      sumpartmatL[, j] <- f1_minus_f2
      omprvecL[j] <- omega_der
      
      b <- b_mat_bon[j, 2]
      C1 <- C[j]
      C2 <- C[J]
      g1 <- gamma[j]
      g2 <- gamma[J]

      f1_minus_f2 <- pos(b - C1 * Norm(Vplus(X, mon_ind))^g1 - 
                           C2 * Norm(Vminus(X, mon_ind))^g2) / sigma      
      omega_der_denom <- sum(f1_minus_f2 / sigma)
      omega_der <- delta / omega_der_denom
      
      sumpartmatU[, j] <- f1_minus_f2
      omprvecU[j] <- omega_der
      
    }
    
    omprvec <- c(omprvecL, omprvecU)
    sumpartmat <- cbind(sumpartmatL, -sumpartmatU)
    
    muvec <- stats::qnorm(a) * omprvec
    sigmamat <- matrix(0, 2 * J, 2 * J)
    
    for (j1 in 1:(2 * J)) {
      for (j2 in j1:(2 * J)) {
        if (j1 == j2) {
          sigmamat[j1, j1] <- omprvec[j1]^2
        } else {
          sigmamat[j1, j2] <- (omprvec[j1] * omprvec[j2] / delta^2) *
            sum(sumpartmat[, j1] * sumpartmat[, j2])
          sigmamat[j2, j1] <- sigmamat[j1, j2]
        }
      }
    }
    
    q <- maxminQ(muvec, sigmamat, alpha, simlen, lower = T)
    return(q)
    
  }
  
  r <- stats::uniroot(f, c(alpha / (2 * J), alpha))
  res <- r$root
  
  return(res) 
}

##' (One-sided) Optimal "tau"
##'
##' Calculates the optimal "tau" 
##' 
##' @param gamma length two vector of gammas (\eqn{(\gamma_1, \gamma_2)'})
##' @param C length two vector of Cs (\eqn{(C_1, C_2)'})
##' @param gam_min
##' @param C_max
##' @param Xt \eqn{n_t \times k} design matrix for the treated units
##' @param Xc \eqn{n_c \times k} design matrix for the control units
##' @param mon_ind indice of the monotone variables
##' @param sigma_t standard deviation of the error term for the treated units
##' (either length 1 or \eqn{n_t})
##' @param sigma_c standard deviation of the error term for the control units
##' @param lower 
##' @param simlen 
##' @param alpha 
##' @param Yt length \eqn{n_t} of outcome variables for the treated units
##' @param Yc length \eqn{n_c} of outcome variables for the control units
##' @export
AdjAlpha_RD <- function(gamma, C, gam_min = min(gamma), C_max = max(C), Xt, Xc, 
                        mon_ind, sigma_t, sigma_c, lower, simlen = 1e05, alpha = .05){
  
  J <- length(gamma)
  nt <- nrow(Xt)
  nc <- nrow(Xc)
  
  prob_maxV <- function(tau){
    
    delta <- stats::qnorm(1 - tau)
    b_mat_bon <- matrix(0, J, 2)
    
    weightmat_t <- matrix(0, nt, J)
    weightmat_c <- matrix(0, nc, J)
    
    for (j in 1:J) {

      gampair_j <- c(gamma[j], gam_min)
      Cpair_j <- c(C[j], C_max)
      swap <- lower
     
      modres <- modsol_RD(delta, gampair_j, Cpair_j, Xt, Xc, mon_ind,
                          sigma_t, sigma_c, swap = swap)
      
      bt <- modres$bt
      bc <- modres$bc

      K_tj <- K_fun(bt, gamma = gampair_j, C = Cpair_j, X = Xt,
                    mon_ind = mon_ind, swap = swap)
      K_cj <- K_fun(bc, gamma = gampair_j, C = Cpair_j, X = Xc,
                    mon_ind = mon_ind, swap = !swap)
      
      weightmat_t[, j]  <- bt * K_tj / sigma_t / delta
      weightmat_c[, j]  <- bc * K_cj / sigma_c / delta  
    }
   
    sigmamat_t <- t(weightmat_t) %*% weightmat_t
    sigmamat_c <- t(weightmat_c) %*% weightmat_c
    sigmamat <- sigmamat_t + sigmamat_c

    q <- maxminQ2(numeric(J), sigmamat, alpha, simlen, tau)

    return(q) 
  }
  
  r <- stats::uniroot(prob_maxV, c(alpha / J, alpha), extendInt = "yes")
  res <- r$root
  
  return(res)
}
