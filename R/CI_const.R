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

  Kt <- K_fun(bt, rep(gamma,2), rep(C,2), Xt, mon_ind)
  gf_ip_iota_t <- sum(bt * Kt / sigma_t^2)

  Kc <- K_fun(bc, rep(gamma,2), rep(C,2), Xc, mon_ind)
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

##' Minimax CI for the RD parameter
##'
##' This function calculates shortest minimax (fixed-length) CI for the RD
##' parameter
##' 
##' @param Yt 
##' @param Yc 
##' @param Xt 
##' @param Xc 
##' @param gam_min 
##' @param C_max 
##' @param mon_ind 
##' @param sigma_t 
##' @param sigma_c 
##' @param alpha 
##' @param opt_b 
##' @param min_half_length 
##' @param maxb.const 
##' @param Prov.Plot 
##'
##' @export
CI_minimax_RD <- function(Yt, Yc, Xt, Xc, gam_min, C_max, mon_ind, sigma_t, sigma_c,
                          alpha = .05, opt_b = NULL, min_half_length = NULL,
                          maxb.const = 10, Prov.Plot = FALSE) {
  
  if(is.null(opt_b) | is.null(min_half_length)){
    
    modres <- modsol_RD(0,rep(gam_min,2), rep(C_max,2), Xt, Xc, mon_ind, 
                        sigma_t, sigma_c)
    minbt <- modres$bt
    minbc <- modres$bc
    minb <- minbt + minbc
    
    maxb <- maxb.const * (modsol_RD(qnorm(1 - alpha/2),rep(gam_min,2), rep(C_max,2), 
                                   Xt, Xc, mon_ind, sigma_t, sigma_c)$bt +
                            modsol_RD(qnorm(1 - alpha/2),rep(gam_min,2), rep(C_max,2), 
                                      Xt, Xc, mon_ind, sigma_t, sigma_c)$bc)
    
    CI_length_sol <- stats::optimize(CI_length_RD, interval = c(minb, maxb), 
                              gamma = gam_min, C = C_max,
                              Xt = Xt, Xc = Xc, mon_ind = mon_ind,
                              sigma_t = sigma_t, sigma_c = sigma_c, alpha = alpha)
    
    min_half_length <- CI_length_sol$objective / 2
    opt_b <- CI_length_sol$minimum
    
    if(Prov.Plot == TRUE){
      
      numgrid = 100
      xintv = seq(from = minb, to = maxb, length.out = numgrid)
      yvec = numeric(numgrid)
      for(i in 1:numgrid){
        
        yvec[i] = CI_length_RD(xintv[i], gamma = gam_min, C = C_max,
                              Xt = Xt, Xc = Xc, mon_ind = mon_ind,
                              sigma_t = sigma_t, sigma_c = sigma_c, alpha = alpha)
      }
      
      plot(xintv, yvec, type = "l", xlab = "modulus", ylab = "CI_length")
      abline(v = opt_b, col = "red", lty = 2)
    }
  }
  
  opt_Lhat <- Lhat_RD_fun(b = opt_b, gamma = rep(gam_min, 2), C = rep(C_max, 2),
                          Xt = Xt, Xc = Xc, Yt = Yt, Yc = Yc, mon_ind = mon_ind,
                          sigma_t = sigma_t, sigma_c = sigma_c)
  res = c(opt_Lhat - min_half_length, opt_Lhat + min_half_length)

  return(res)
}


##' Minimax CI Modulus calcuation for the RD parameter
##'
##' This function calculates the length of fixed length minimax CI and 
##' optimal modulus value for the RD parameter
#'
#' @param Xt 
#' @param Xc 
#' @param gam_min 
#' @param C_max 
#' @param mon_ind 
#' @param sigma_t 
#' @param sigma_c 
#' @param alpha 
#' @param maxb.const 
#' @param Prov.Plot 
#' 
#' @export
CI_minimax_RD_mod <- function(Xt, Xc, gam_min, C_max, mon_ind, sigma_t, sigma_c,
                          alpha = .05, maxb.const = 10, Prov.Plot = FALSE) {
  
  modres <- modsol_RD(0,rep(gam_min,2), rep(C_max,2), Xt, Xc, mon_ind, 
                      sigma_t, sigma_c)
  minbt <- modres$bt
  minbc <- modres$bc
  minb <- minbt + minbc
  
  maxb <- maxb.const * (modsol_RD(qnorm(1 - alpha/2),rep(gam_min,2), rep(C_max,2), 
                                  Xt, Xc, mon_ind, sigma_t, sigma_c)$bt +
                          modsol_RD(qnorm(1 - alpha/2),rep(gam_min,2), rep(C_max,2), 
                                    Xt, Xc, mon_ind, sigma_t, sigma_c)$bc)
  
  CI_length_sol <- stats::optimize(CI_length_RD, interval = c(minb, maxb), 
                                   gamma = gam_min, C = C_max,
                                   Xt = Xt, Xc = Xc, mon_ind = mon_ind,
                                   sigma_t = sigma_t, sigma_c = sigma_c, alpha = alpha)
  
  min_half_length <- CI_length_sol$objective / 2
  opt_b <- CI_length_sol$minimum
  
  if(Prov.Plot == TRUE){
    
    numgrid = 100
    xintv = seq(from = minb, to = maxb, length.out = numgrid)
    yvec = numeric(numgrid)
    for(i in 1:numgrid){
      
      yvec[i] = CI_length_RD(xintv[i], gamma = gam_min, C = C_max,
                             Xt = Xt, Xc = Xc, mon_ind = mon_ind,
                             sigma_t = sigma_t, sigma_c = sigma_c, alpha = alpha)
    }
    
    plot(xintv, yvec, type = "l", xlab = "modulus", ylab = "CI_length")
    abline(v = opt_b, col = "red", lty = 2)
  }
  
  res = list(opt_b = opt_b, min_half_length = min_half_length)
  
  return(res)
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
##'
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
##' @param corr_tol
##'
##' @export
AdjAlpha_RD <- function(gamma, C, gam_min = min(gamma), C_max = max(C), Xt, Xc, 
                        mon_ind, sigma_t, sigma_c, lower, simlen = 1e05, alpha, 
                        corr_tol = 1 - 10^(-3)){
  
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
    
    corrs <- sigmamat[row(sigmamat) != col(sigmamat)]
    if(min(corrs) > corr_tol){
      
      q <- 0
    }else{
      
      q <- maxminQ2(numeric(J), sigmamat, alpha, simlen, tau)
    }

    return(q) 
  }
  
  if(prob_maxV(alpha) == 0){
    
    res <- alpha
  }else{
    
    r <- stats::uniroot(prob_maxV, c(alpha / J, alpha), extendInt = "yes")
    res <- r$root
  }
  
  return(res)
}


##' (One-sided) Procedure choice
##'
##' Calculates the adaptive procedure "closest" to the oracle one
##' 
##' @param Cgrid
##' @param gamma length two vector of gammas (\eqn{(\gamma_1, \gamma_2)'})
##' @param C length two vector of Cs (\eqn{(C_1, C_2)'})
##' @param Xt \eqn{n_t \times k} design matrix for the treated units
##' @param Xc \eqn{n_c \times k} design matrix for the control units
##' @param mon_ind indice of the monotone variables
##' @param sigma_t standard deviation of the error term for the treated units
##' (either length 1 or \eqn{n_t})
##' @param sigma_c standard deviation of the error term for the control units
##' @param lower 
##' @param alpha 
##' @param procmat
##' @param Procnames
##' @param Nsim_proc
##' @param rhofun corresponds to \eqn{\rho} in our documentation
##' @export

which_proc <- function(Cgrid, gamma, C, Xt, Xc, mon_ind, sigma_t, sigma_c,
                       lower, alpha, procmat, ProcNames, Nsim_proc, 
                       rhofun = c("diff","ratio")){
  
  rhofun <- match.arg(rhofun)
  
  if(rhofun == "diff"){
    
    rhofun = function(l1, l2){  
      
      res <- l1 - l2
      
      return(res)
    }
  }else{
    
    rhofun = function(l1, l2){  
      
      res <- l1 / l2
      
      return(res)
    }
  }
  
  gamgrid <- rep(1, length(Cgrid))
  num_grid <- length(Cgrid)
  N_proc <- length(ProcNames)
  
  gam_min <- min(gamgrid)
  C_max <- max(Cgrid)

  nt <- nrow(Xt)   # Number of treated sample
  nc <- nrow(Xc)  # Number of control sample
  k <- ncol(Xt)
  
  fXt_wc <- matrix(NA, nrow = nt, ncol = num_grid)
  fXc_wc <- matrix(NA, nrow = nc, ncol = num_grid)
  
  for(j in 1:num_grid) {

    fXt_wc[, j] <- - Cgrid[j] * Norm(Vminus(Xt, mon_ind))^gamgrid[j] 
    fXc_wc[, j] <- Cgrid[j] * Norm(Vplus(Xc, mon_ind))^gamgrid[j]
  }
  
  # Calculating moduli of continuity and related quantities
  
  proc_bmat = matrix(0, sum(procmat), 4)
  proc_dmat = matrix(0, sum(procmat), 4)
  proc_alphavec = numeric(N_proc)
  
  beg_ind <- 1
  
  for(i in 1:N_proc){
    
    print(paste("proc",i))
    
    adpt_ind <- procmat[i, ]
    adpt_len <- sum(procmat[i, ])
    end_ind <- beg_ind + adpt_len -1
    modres <- mod_del_cal(gamma[adpt_ind], C[adpt_ind], gam_min, C_max, 
                          Xt, Xc, mon_ind, sigma_t, sigma_c)
    proc_bmat[beg_ind:end_ind, ] <- modres$b_mat
    proc_dmat[beg_ind:(beg_ind + adpt_len -1), ] <- modres$delta_mat
    
    if(lower == T){
      proc_alphavec[i] <- modres$alnewL
    }else{
      proc_alphavec[i] <- modres$alnewU
    }
    
    
    beg_ind <- beg_ind + adpt_len
    
  }
  
  CIlen_orc <- numeric(num_grid)
  
  for(j in 1:num_grid){
    # res_mod_orc <- mod_del_cal_orc(c(gamgrid[j], min(gamgrid)), c(Cgrid[j], max(Cgrid)),
    #                                Xt, Xc, mon_ind, sigma_t, sigma_c)
    # CIlen_orc[j] <- ifelse(lower, sum(res_mod_orc$b_vec[1:2]), sum(res_mod_orc$b_vec[3:4]))
    res_mod_orc <- mod_del_cal(gamgrid[j], Cgrid[j], gam_min, C_max,
                                   Xt, Xc, mon_ind, sigma_t, sigma_c)
    CIlen_orc[j] <- ifelse(lower, sum(res_mod_orc$b_mat[, 1:2]), 
                           sum(res_mod_orc$b_mat[, 3:4]))
  }
  
  # Generating y_i's and computating CIs
  
  CIs <- array(0, c(num_grid, Nsim_proc, N_proc))
  # CIs_orc <- matrix(0, Nsim_proc, num_grid)
  
  
  for (i in 1:Nsim_proc) {
    
    ut <- stats::rnorm(nt, sd = sigma_t)   # generate u_i's for the treated
    Yt <- fXt_wc + ut

    uc <- stats::rnorm(nc, sd = sigma_c)   # generate u_i's for the contol
    Yc <- fXc_wc + uc
    
    for (j in 1:num_grid){
      
      beg_ind <- 1
      
      for (j2 in 1:N_proc){
        
        adpt_ind <- procmat[j2,]
        adpt_len <- sum(procmat[j2,])
        end_ind <- beg_ind + adpt_len -1
        
        CIs[j,i,j2] <- CI_one_sd_RD(proc_bmat[beg_ind:end_ind, , drop=F],
                                    proc_dmat[beg_ind:end_ind, , drop=F],
                                    gamma[adpt_ind], C[adpt_ind], gam_min, C_max,
                                    Xt, Xc, mon_ind, sigma_t, sigma_c, Yt[,j], Yc[,j],
                                    lower, alpha = proc_alphavec[j2])
        
        beg_ind <- beg_ind + adpt_len
      }
    }
    
    if(i%%50 == 0) print(i)
  }
  
  CIlen <- apply(0 - CIs,MARGIN = c(1,3),FUN=mean,na.rm=T)
  
  if(lower == T){
    
    CIcov <- apply(CIs < 0,MARGIN = c(1,3),FUN=mean,na.rm=T)
  }else{
    
    CIcov <- apply(CIs > 0,MARGIN = c(1,3),FUN=mean,na.rm=T)
  }
  
  if(lower == F){
    
    CIlen <- -CIlen
  }
  
  regret_mat <- rhofun(CIlen, 
                       matrix(rep(CIlen_orc, N_proc), num_grid, N_proc, byrow = F))
  
  maxregret <- apply(regret_mat, 2, max)
  names(maxregret) <- ProcNames
  
  minMRind <- which.min(maxregret)
  bestProc <- C[procmat[minMRind, ]]
  bestProcName <- ProcNames[minMRind]
  
  print(paste("C", seq(1:length(bestProc)), ":", bestProc), sep="")
  
  res <- list(CIlen = CIlen, CIlen_orc = CIlen_orc, regret_mat = regret_mat,
              maxregret = maxregret, bestProc = bestProc, bestProcName = bestProcName,
              proc_bmat = proc_bmat, proc_dmat = proc_dmat,
              proc_alphavec = proc_alphavec, CIcov = CIcov)
  
  return(res)
}


#' Title
#'
#' @param Yt 
#' @param Yc 
#' @param Xt 
#' @param Xc 
#' @param C 
#' @param C_max 
#' @param mon_ind 
#' @param sigma_t 
#' @param sigma_c 
#' @param lower 
#' @param two_sided 
#' @param alpha 
#' @param modres 
#'
#' @return
#' @export
#'
#' @examples
CI_gen <- function(Yt, Yc, Xt, Xc, C, C_max = max(C), mon_ind, sigma_t, sigma_c,
                   lower, two_sided = F, alpha,
                   modres = NULL){

  gamma <- rep(1, length(C))

  alpha_new <- ifelse(two_sided, alpha / 2, alpha)

  if(is.null(modres)){

    modres <- mod_del_cal(gamma, C, 1, C_max, Xt, Xc, mon_ind, sigma_t, sigma_c, 
                          alpha_new)
  }

  bmat <- modres$b_mat
  dmat <- modres$delta_mat

  if(lower == T){

    al_adj <- modres$alnewL
  }else{

    al_adj <- modres$alnewU
  }

  res1 <- CI_one_sd_RD(bmat, dmat, gamma, C, 1, C_max, Xt, Xc, mon_ind, sigma_t,
                       sigma_c, Yt, Yc, lower, alpha = al_adj)
  
  if(two_sided == F){
    
    res_l <- ifelse(lower, res1, -Inf)
    res_u <- ifelse(lower, Inf, res1)
    
  }else{

    if(length(C) > 1){

      modres_mm <- mod_del_cal(1, C_max, 1, C_max, Xt, Xc, mon_ind,
                               sigma_t, sigma_c, alpha_new)
      
    }else{

      modres_mm <- modres
    }
    
    bmat_mm <- modres_mm$b_mat
    deltamat_mm <- modres_mm$delta_mat
    
    res2 <- CI_one_sd_RD(bmat_mm, deltamat_mm, 1, C_max, 1, C_max, Xt, Xc, 
                         mon_ind, sigma_t, sigma_c, Yt, Yc, lower = !lower, 
                         alpha_new)
    
    res_l <- ifelse(lower, res1, res2)
    res_u <- ifelse(lower, res2, res1)

  }

  res <- list(CI_l = res_l, CI_u = res_u, alpha = alpha, Cvec = C, C_max = C_max)
  
  return(res)
}