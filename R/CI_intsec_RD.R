cov_calc <- function(delta, ht_j, hc_j, ht_k, hc_k, Cj, Ck, Cbar, Xt, Xc, mon_ind, 
                       sigma_t, sigma_c){
  
  if(missing(ht_j) | missing(hc_j) | missing(ht_k) | missing(hc_k)){
    
    hres_j <- bw_adpt(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
    ht_j <- hres_j$ht
    hc_j <- hres_j$hc
    
    hres_k <- bw_adpt(delta, Ck, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
    ht_k <- hres_k$ht
    hc_k <- hres_k$hc
  }
  
  num_tj <- K_fun(b = ht_j, C_pair = c(Cbar, Cj), X = Xt, mon_ind = mon_ind)
  num_tk <- K_fun(b = ht_k, C_pair = c(Cbar, Ck), X = Xt, mon_ind = mon_ind)
  num_t <- sum(num_tj * num_tk / sigma_t^2)
  sum_t <- ht_j * ht_k * num_t / delta^2
  
  num_cj <- K_fun(b = hc_j, C_pair = c(Cj, Cbar), X = Xc, mon_ind = mon_ind)
  num_ck <- K_fun(b = hc_k, C_pair = c(Ck, Cbar), X = Xc, mon_ind = mon_ind)
  num_c <- sum(num_cj * num_ck / sigma_c^2)
  sum_c <- hc_j * hc_k * num_c / delta^2
  
  res <- sum_t + sum_c
  return(res)
}

cov_mat_calc <- function(delta, Cvec, Xt, Xc, mon_ind, sigma_t, sigma_c, hmat){
  
  J <- length(Cvec)
  Cbar <- max(Cvec)
  
  if(missing(hmat)){
    
    hmat <- matrix(0, nrow = J, ncol = 2)
    
    for(j in 1:J){
      
      Cj <- Cvec[j]
      
      hres_j <- bw_adpt(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
      hmat[j, 1] <- hres_j$ht
      hmat[j, 2] <- hres_j$hc
    }
  }
  
  res <- diag(1, J, J)
  
  for(j in 1:(J - 1)){
    
    for(k in (j + 1):J){
      
      Cj <- Cvec[j]
      Ck <- Cvec[k]
      ht_j <- hmat[j, 1]
      hc_j <- hmat[j, 2]
      ht_k <- hmat[k, 1]
      hc_k <- hmat[k, 2]
      
      res[j, k] <- cov_calc(delta, ht_j, hc_j, ht_k, hc_k, Cj, Ck, Cbar, 
                             Xt, Xc, mon_ind, sigma_t, sigma_c)
      res[k, j] <- res[j, k]
    }
  }
  
  return(res)
}

max_Q <- function(covmat, alpha, num_sim = 10^4){
  
  J <- nrow(covmat)
  rn <- MASS::mvrnorm(n = num_sim, mu = rep(0, J), Sigma = covmat)
  rn_max <- apply(rn, 1, max)
  
  res <- quantile(rn_max, 1 - alpha)
  return(res)
}

tau_star <- function(Cvec, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha, num_sim = 10^4,
                     grid_len = 100){
  
  eqn <- function(tau){
    
    delta <- qnorm(1 - tau)
    covmat <- cov_mat_calc(delta, Cvec, Xt, Xc, mon_ind, sigma_t, sigma_c)
    
    res <- max_Q(covmat, alpha, num_sim) - qnorm(1 - tau)
    return(res)
  }
  
  J <- length(Cvec)
  grid <- seq(from = alpha / J, to = alpha * 1.01, length.out = grid_len)
  eqn_res <- numeric(grid_len) 
  
  for(i in 1:grid_len){
    
    eqn_res[i] <- eqn(grid[i])
  }
  
  eqn_res <- abs(sort(eqn_res))
  tau_ind <- which.min(eqn_res)
  
  res <- grid[tau_ind]
  return(res)
}
