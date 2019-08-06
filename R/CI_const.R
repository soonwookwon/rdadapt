CI_length <- function(b, gamma, C, X, mon_ind, sigma = 1) {
  
  om_inv <- invmod(b, rep(gamma, 2), rep(C, 2), X, mon_ind, sigma)

  if (om_inv == 0) return(Inf)
  
  gf_ip_iota <- pos(b - C * Norm(Vplus(X,mon_ind))^gamma -
                      C * Norm(Vminus(X,mon_ind))^gamma) / sigma
  
  gf_ip_iota <- sum(gf_ip_iota)
  
  sd <- om_inv / gf_ip_iota
  bias <- .5 * ( b - ( om_inv^2 / gf_ip_iota ))
  
  cva <- ifelse((bias / sd) > 3,
                (bias / sd) + qnorm(1 - alpha),
                sqrt(qchisq(1 - alpha, df = 1, ncp = (bias / sd)^2)))
  
  return( 2 * cva * sd)
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
                      sigma_t = sigma_t, sigma_c = sigma_c, ret_list = TRUE)

  bc <- om_inv$bc
  bt <- om_inv$bt
  ominv_t <- om_inv$dt 
  ominv_c <- om_inv$dc
  om_inv <- sqrt((ominv_t)^2 + (ominv_c)^2)
    
  if (om_inv == 0) return(Inf)
  
  gf_ip_iota_t <- pos(bt - C * Norm(Vplus(Xt,mon_ind))^gam -
                        C * Norm(Vminus(Xt,mon_ind))^gam) / sigma_t
  gf_ip_iota_t <- sum(gf_ip_iota_t)
  
  gf_ip_iota_c <- pos(bc - C * Norm(Vplus(Xc,mon_ind))^gam -
                        C * Norm(Vminus(Xc,mon_ind))^gam) / sigma_c
  gf_ip_iota_c <- sum(gf_ip_iota_c)
  
  sd <- sqrt((ominv_t/gf_ip_iota_t)^2 + (ominv_c/gf_ip_iota_c)^2)
  
  bias <- .5 * (b - (om_inv^2 /  (gf_ip_iota_t + gf_ip_iota_c)))
  cva <- ifelse(abs(bias / sd) > 3,
                abs(bias / sd) + qnorm(1 - alpha),
                sqrt(qchisq(1 - alpha, df = 1, ncp = (bias / sd)^2)))
  
  return(2 * cva * sd)
}


minimax_Lhat <- function(b, gamma, C, X, Y, mon_ind, sigma) {
  
  f1 <- C * Norm(Vplus(X, mon_ind))^gamma
  f2 <- pmax(b - C * Norm(Vminus(X, mon_ind))^gamma, f1)
    
  return(.5 * b + sum((f2 - f1) * (Y - .5 * (f1 + f2)) / sigma^2) /
           sum((f2 - f1) / sigma^2))
}

minimax_Lhat_RD <- function(b, gamma, C, Xt, Xc, Yt, Yc, sigma_t, sigma_c,
                            mon_ind) {
  
  om_inv <- invmod_RD(b, rep(gamma, 2), rep(C, 2), Xt, Xc, mon_ind = mon_ind,
                      sigma_t = sigma_t, sigma_c = sigma_c, ret_list = TRUE)
  
  bc <- om_inv$bc
  bt <- om_inv$bt
  
  ft1 <- C * Norm(Vplus(Xt, mon_ind))^gamma
  ft2 <- pmax(bt - C * Norm(Vminus(Xt,mon_ind))^gamma, ft1)
  
  fc1 <- C * Norm(Vplus(Xc, mon_ind))^gamma
  fc2 <- pmax(bc - C * Norm(Vminus(Xc, mon_ind))^gamma, fc1)

  Lhat_t <- .5 * bt + sum((ft2 - ft1) * (Yt - .5 * (ft1 + ft2)) / sigma_t^2) /
    sum((ft2 - ft1) / sigma_t^2)

  Lhat_c <- .5 * bc + sum((fc2 - fc1) * (Yc - .5 * (fc1 + fc2)) / sigma_c^2) /
    sum((fc2 - fc1) / sigma_c^2)

  return(Lhat_t - Lhat_c)
}


CI_one_sd <- function(bpairmat, delta, gamma, C, X, mon_ind, sigma, y, lower, al,
                      Dir = FALSE, maxQ = FALSE, simlen = 5e04, qres = FALSE){
  
  # bpairmat[j,] = (omega(del,F_J,F_j, omega(del,F_j,F_J)) 
  # where delta is given by sigma*(z_beta + z_{1-al/(2J)})
  
  J <- length(gamma)
  n <- length(y)
  
  if (Dir) {
    endj <- 1
    chatvec <- numeric(1)
  } else {
    endj <- J
    chatvec <- numeric(J)
  }
  
  if (maxQ) {
    sumpartmat <- matrix(0, n, J)
    omprvec <- numeric(J)
  } 
  
  for (j in 1:endj) {
    
    if (lower == TRUE) {
      b <- bpairmat[j, 1]
      C1 <- C[J]
      C2 <- C[j]
      g1 <- gamm[J]
      g2 <- gamma[j]
    } else {
      b <- bpairmat[j, 2]
      C1 <- C[j]
      C2 <- C[J]
      g1 <- gamma[j]
      g2 <- gamma[J]
    } 
    
    f1_minus_f2 <- pos(b - C1 * Norm(Vplus(X, mon_ind))^g1 - 
                         C2 * Norm(Vminus(X, mon_ind))^g2) / sigma
    y_minus_midpoint <- (y - (C1 * Norm(Vplus(X,mon_ind))^g1 + b - 
                                C2 * Norm(Vminus(X,mon_ind))^g2) / 2) / sigma
    
    Lhat_frac_num <- sum(f1_minus_f2 * y_minus_midpoint)
    
    # Calculation of omega'(del)
    # Uses Lemma B.3 in A&K(2016)

    omega_der_denom <- sum(f1_minus_f2 / sigma)
    
    omega_der <- delta / omega_der_denom 
    
    Lhat <- .5 * b + Lhat_frac_num / omega_der_denom
    
    Ubias <- .5 * (b - delta * omega_der)
    # chat <- ifelse(lower == TRUE,
    #               Lhat - Ubias - qnorm(1-al)*sigma*ompr,
    #               -(Lhat - Ubias - qnorm(1-al)*sigma*ompr))
    
    chat <- ifelse(lower == TRUE,
                  Lhat - Ubias - qnorm(1 - al) * omega_der,
                  Lhat + Ubias + qnorm(1 - al) * omega_der)
    
    chatvec[j] <- chat
    
    if (maxQ == T) {
      sumpartmat[, j] <- f1_minus_f2
      omprvec[j] <- omega_der
    } 
    
  }
  
  if (maxQ == T) {
    
    muvec <- ifelse(lower==T, qnorm(al), qnorm(1 - al)) * omprvec
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
    
    q <- maxminQ(muvec, sigmamat, alpha / 2, simlen, lower)
    res <- ifelse(lower == T, max(chatvec) - q, min(chatvec) - q)

  } else {
    res <- ifelse(lower == T,max(chatvec),min(chatvec))
  }
  
  if (qres == T) {
    res <- c(res,q)
  }
  
  return(res)
}


CI_one_sd_RD <- function(bpairmat, dmat, gamma, C, Xt, Xc, mon_ind,
                         sigma_t, sigma_c, yt, yc, lower, al, Dir = F){
  
  # bpairmat[j,] = (omega(del,F_J,F_j, omega(del,F_j,F_J)) 
  # where del is given by sigma*(z_beta + z_{1-al/(2J)})
  
  J <- length(gam)
  
  bpairmat_l <- bpairmat[,1:2,drop=F]
  bpairmat_u <- bpairmat[,3:4,drop=F]
  dmat_l <- dmat[,1:2,drop=F]
  dmat_u <- dmat[,3:4,drop=F]
  
  if(Dir == TRUE){
    endj <- 1
    chatvec <- numeric(1)
  }else{
    endj <- J
    chatvec <- numeric(J)
  }
  
  for(j in 1:endj){
    
    if(lower == TRUE){
      bt <- bpairmat_l[j,1]
      bc <- bpairmat_l[j,2]
      del_t <- dmat_l[j,1]
      del_c <- dmat_l[j,2]
      C1 <- C[J]
      C2 <- C[j]
      g1 <- gam[J]
      g2 <- gam[j]
    } 
    else{
      bt <- bpairmat_u[j,1]
      bc <- bpairmat_u[j,2]
      del_t <- dmat_u[j,1]
      del_c <- dmat_u[j,2]
      C1 <- C[j]
      C2 <- C[J]
      g1 <- gam[j]
      g2 <- gam[J]
    } 
    
    Lhatfun <- function(b,C1,C2,g1,g2,X,y){
      sumpart1 <- (b - C1 * Norm(Vplus(X,mon_ind))^g1 - 
                    C2 * Norm(Vminus(X,mon_ind))^g2) *  
        (b - C1 * Norm(Vplus(X,mon_ind))^g1 - 
           C2 * Norm(Vminus(X,mon_ind))^g2 > 0) 
      sumpart2 <- y - (C1 * Norm(Vplus(X,mon_ind))^g1 + b - 
                        C2 * Norm(Vminus(X,mon_ind))^g2)/2
      
      sumpart <- sumpart1 * sumpart2
      
      denom <- sum((b - C1 * Norm(Vplus(X,mon_ind))^g1 - 
                     C2 * Norm(Vminus(X,mon_ind))^g2) *  
                    (b - C1 * Norm(Vplus(X,mon_ind))^g1 - 
                       C2 * Norm(Vminus(X,mon_ind))^g2 > 0)) 
      
      Lhat <- b/2 + sum(sumpart) / denom
      
      res <- c(Lhat,denom) # denom is used to compute omega'
      
      return(res)
    }
    
    Lhat_t_res <- Lhatfun(bt,C1,C2,g1,g2,Xt,yt)
    Lhat_c_res <- Lhatfun(bc,C2,C1,g2,g1,Xc,yc)
    
    ompr_t <- del_t / Lhat_t_res[2]
    ompr_c <- del_c / Lhat_c_res[2]
    ompr <- sqrt(ompr_t^2 + ompr_c^2)
    
    Lhat <- Lhat_t_res[1] - Lhat_c_res[1]
    
    chat <- ifelse(lower == TRUE,
                  Lhat - 0.5*((bt + bc) + qnorm(1-al)*ompr),
                  Lhat + 0.5*((bt + bc) + qnorm(1-al)*ompr))
    
    
    chatvec[j] <- chat
    
    
  }
  
  res <- ifelse(lower == T,max(chatvec),min(chatvec))
  
  return(res)
  
}


AdjAlpha <- function(gam,C,X,mon_ind,lower,simlen = 10000){
  
  J <- length(gam)
  n <- length(X[,1])
  
  f <- function(a){
    
    del_adpt <- qnorm(1-a)
    b_mat_bon <- matrix(0,J,2)
    
    for(j in 1:J){
      gampair_j <- c(gam[j],gam[J])
      Cpair_j <- c(C[j],C[J])
      
      b_mat_bon[j,1] <- modsol(del_adpt,gampair_j[2:1],Cpair_j[2:1],X,mon_ind)
      b_mat_bon[j,2] <- modsol(del_adpt,gampair_j,Cpair_j,X,mon_ind)
    }
    
    sumpartmat <- matrix(0,n,J)
    omprvec <- numeric(J)
    
    del <- del_adpt
    
    for(j in 1:J){
      
      if(lower == TRUE){
        b <- b_mat_bon[j,1]
        C1 <- C[J]
        C2 <- C[j]
        g1 <- gam[J]
        g2 <- gam[j]
      } 
      else{
        b <- b_mat_bon[j,2]
        C1 <- C[j]
        C2 <- C[J]
        g1 <- gam[j]
        g2 <- gam[J]
      } 
      
      sumpart1 <- (b - C1 * Norm(Vplus(X,mon_ind))^g1 - 
                    C2 * Norm(Vminus(X,mon_ind))^g2) *  
        (b - C1 * Norm(Vplus(X,mon_ind))^g1 - 
           C2 * Norm(Vminus(X,mon_ind))^g2 > 0) 
      
      ompr <- del / sum(sumpart1)
      sumpartmat[,j] <- sumpart1
      omprvec[j] <- ompr
      
    }
      
    muvec <- ifelse(lower==T, qnorm(a), qnorm(1-a)) * omprvec
    sigmamat <- matrix(0,J,J)
    
    for(j1 in 1:J){
      for(j2 in j1:J){
        if(j1 == j2){
          sigmamat[j1,j1] <- omprvec[j1]^2
        }else{
          sigmamat[j1,j2] <- (omprvec[j1] * omprvec[j2] / del^2) *
            sum(sumpartmat[,j1] * sumpartmat[,j2])
          sigmamat[j2,j1] <- sigmamat[j1,j2]
        }
      }
    }
      
    q <- maxminQ(muvec,sigmamat,alpha/2,simlen,lower)
    return(q)
    
  }
    
  r <- uniroot(f,c(alpha/(2*J),alpha/2))
  res <- r$root
  
  return(res)
  
}



# Adjust alpha for lower and upper together
AdjAlpha2 <- function(gam,C,X,mon_ind,simlen = 10000){
  
  J <- length(gam)
  n <- length(X[,1])
  
  f <- function(a){
    
    del_adpt <- qnorm(1-a)
    b_mat_bon <- matrix(0,J,2)
    
    for(j in 1:J){
      gampair_j <- c(gam[j],gam[J])
      Cpair_j <- c(C[j],C[J])
      
      b_mat_bon[j,1] <- modsol(del_adpt,gampair_j[2:1],Cpair_j[2:1],X,mon_ind)
      b_mat_bon[j,2] <- modsol(del_adpt,gampair_j,Cpair_j,X,mon_ind)
    }
    
    sumpartmatL <- matrix(0,n,J)
    omprvecL <- numeric(J)
    sumpartmatU <- matrix(0,n,J)
    omprvecU <- numeric(J)
    
    del <- del_adpt
    
    for(j in 1:J){
      
      b <- b_mat_bon[j,1]
      C1 <- C[J]
      C2 <- C[j]
      g1 <- gam[J]
      g2 <- gam[j]
      
      sumpart1 <- (b - C1 * Norm(Vplus(X,mon_ind))^g1 - 
                    C2 * Norm(Vminus(X,mon_ind))^g2) *  
        (b - C1 * Norm(Vplus(X,mon_ind))^g1 - 
           C2 * Norm(Vminus(X,mon_ind))^g2 > 0) 
      
      ompr <- del / sum(sumpart1)
      sumpartmatL[,j] <- sumpart1
      omprvecL[j] <- ompr
      
      b <- b_mat_bon[j,2]
      C1 <- C[j]
      C2 <- C[J]
      g1 <- gam[j]
      g2 <- gam[J]
      
      sumpart1 <- (b - C1 * Norm(Vplus(X,mon_ind))^g1 - 
                    C2 * Norm(Vminus(X,mon_ind))^g2) *  
        (b - C1 * Norm(Vplus(X,mon_ind))^g1 - 
           C2 * Norm(Vminus(X,mon_ind))^g2 > 0) 
      
      ompr <- del / sum(sumpart1)
      sumpartmatU[,j] <- sumpart1
      omprvecU[j] <- ompr
      
      
      
    }
    
    omprvec <- c(omprvecL,omprvecU)
    sumpartmat <- cbind(sumpartmatL,-sumpartmatU)
    
    muvec <- qnorm(a) * omprvec
    sigmamat <- matrix(0,2*J,2*J)
    
    for(j1 in 1:(2*J)){
      for(j2 in j1:(2*J)){
        if(j1 == j2){
          sigmamat[j1,j1] <- omprvec[j1]^2
        }else{
          sigmamat[j1,j2] <- (omprvec[j1] * omprvec[j2] / del^2) *
            sum(sumpartmat[,j1] * sumpartmat[,j2])
          sigmamat[j2,j1] <- sigmamat[j1,j2]
        }
      }
    }
    
    q <- maxminQ(muvec,sigmamat,alpha,simlen,Lower=T)
    return(q)
    
  }
  
  r <- uniroot(f,c(alpha/(2*J),alpha))
  res <- r$root
  
  return(res)
  
}


AdjAlpha_RD <- function(gam,C,Xt,Xc,mon_ind,sigma_t,sigma_c,lower,simlen = 10000,rp = alpha){
  
  J <- length(gam)
  nt <- length(Xt[,1])
  nc <- length(Xc[,1])
  
  f <- function(tau){
    
    del <- qnorm(1-tau)
    b_mat_bon <- matrix(0,J,2)
    
    sumpartmat_t <- matrix(0,nt,J)
    sumpartmat_c <- matrix(0,nc,J)
    
    for(j in 1:J){
      
      if(lower==T){
        gampair_j <- c(gam[J],gam[j])
        Cpair_j <- c(C[J],C[j])
      }else if(lower==F){
        gampair_j <- c(gam[j],gam[J])
        Cpair_j <- c(C[j],C[J])
      }
      
      modres <- modsol_RD(del,gampair_j,Cpair_j,Xt,Xc,mon_ind,sigma_t,sigma_c,
                         sol_list = 1)
      
      bt <- modres$bt
      bc <- modres$bc
      C1 <- Cpair_j[1]
      C2 <- Cpair_j[2]
      g1 <- gampair_j[1]
      g2 <- gampair_j[2]
      
      sumpartmat_t[,j]  <- (bt - C1 * Norm(Vplus(Xt,mon_ind))^g1 - 
                    C2 * Norm(Vminus(Xt,mon_ind))^g2) *  
        (bt - C1 * Norm(Vplus(Xt,mon_ind))^g1 - 
           C2 * Norm(Vminus(Xt,mon_ind))^g2 > 0) / sigma_t^2
      
      sumpartmat_c[,j]  <- (bc - C2 * Norm(Vplus(Xc,mon_ind))^g2 - 
                     C1 * Norm(Vminus(Xc,mon_ind))^g1) *  
        (bc - C2 * Norm(Vplus(Xc,mon_ind))^g2 - 
           C1 * Norm(Vminus(Xc,mon_ind))^g1 > 0) / sigma_c^2
      
    }
    

    
    muvec <- rep(0,J)
    sigmamat <- matrix(0,J,J)
    
    for(j1 in 1:J){
      for(j2 in j1:J){
        if(j1 == j2){
          sigmamat[j1,j1] <- 1
        }else{
          sigmamat[j1,j2] <- (sum(sumpartmat_t[,j1] * sumpartmat_t[,j2]) + 
            sum(sumpartmat_c[,j1] * sumpartmat_c[,j2])) / del^2
            
          sigmamat[j2,j1] <- sigmamat[j1,j2]
        }
      }
    }
    
    q <- maxminQ2(muvec,sigmamat,rp,simlen,tau)
    return(q)
    
  }
  
  r <- uniroot(f,c(rp/J,rp))
  res <- r$root
  
  return(res)
  
}
