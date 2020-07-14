#' Calculates the bandwidths h_jt and h_jc
#'
#' @param delta a nonnegative scalar
#' @param Cj the smoothness parameter aiming to adapt to
#' @param Cbar the largest smoothness parameter
#' @param Xt data for the treated
#' @param Xc data for the control
#' @param mon_ind index for monotone variables
#' @param sigma_t standard deviation for the treated group observations
#' @param sigma_c standard deviation for the control group observations
#'
#' @return a list of two values, h_jt and h_jc
#' @export
#'
#' @examples
bw_adpt <- function(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c){
  
  C_pair = c(Cbar, Cj)
  
  modres <- modsol_RD(delta, gam_pair = c(1, 1), C_pair, Xt, Xc, mon_ind, 
                      sigma_t, sigma_c, swap = FALSE)
  
  ht <- modres$bt
  hc <- modres$bc
  res <- list(ht = ht, hc = hc)
  
  return(res)
}

#' Estimator for adaptive CI
#' 
#' Calculates \eqn{\hat{L}_j(\delta)} in the paper
#'
#' @param delta a nonegative scalar value: 
#' it can be left unspecified if `ht` and `hc` are specified
#' @param ht the modulus value for the treated observations: 
#' it can be left unspecified if `delta` is specified
#' @param hc the modulus value for the control observations: 
#' it can be left unspecified if `delta` is specified
#' @param Cj the smoothness parameter aiming to adapt to
#' @param Cbar the largest smoothness parameter
#' @param Xt data for the treated
#' @param Xc data for the control
#' @param mon_ind index for monotone variables
#' @param sigma_t standard deviation for the treated group observations
#' @param sigma_c standard deviation for the control group observations
#' @param Yt outcome value for the treated group observations
#' @param Yc outcome value for the control group observations
#'
#' @return a scalar value of the estimator
#' @export
#'
#' @examples
Lhat_fun_RD <- function(delta, ht, hc, Cj, Cbar, Xt, Xc, mon_ind, 
                     sigma_t, sigma_c, Yt, Yc) {
  
  if(missing(ht) | missing(hc)){
    
    hres <- bw_adpt(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
    ht <- hres$ht
    hc <- hres$hc
  }
  
  num_it <- K_fun(b = ht, gam_pair = c(1, 1), C_pair = c(Cbar, Cj), Xt, mon_ind) * Yt / sigma_t^2
  denom_it <- K_fun(b = ht, gam_pair = c(1, 1), C_pair = c(Cbar, Cj), Xt, mon_ind) / sigma_t^2
  res_t <- sum(num_it) / sum(denom_it)
  
  num_ic <- K_fun(b = hc, gam_pair = c(1, 1), C_pair = c(Cj, Cbar), Xc, mon_ind) * Yc / sigma_c^2
  denom_ic <- K_fun(b = hc, gam_pair = c(1, 1), C_pair = c(Cj, Cbar), Xc, mon_ind) / sigma_c^2
  res_c <- sum(num_ic) / sum(denom_ic)
  
  return(res_t - res_c)
}

#' Worst-case bias calculation helper function
#' 
#' Calculates \eqn{0.5 \times (a_{jt} - a_{jc})} in the worst-case formula
#'
#' @param delta a nonegative scalar value: 
#' it can be left unspecified if `ht` and `hc` are specified
#' @param ht the modulus value for the treated observations: 
#' it can be left unspecified if `delta` is specified
#' @param hc the modulus value for the control observations: 
#' it can be left unspecified if `delta` is specified
#' @param Cj the smoothness parameter aiming to adapt to
#' @param Cbar the largest smoothness parameter
#' @param Xt data for the treated
#' @param Xc data for the control
#' @param mon_ind index for monotone variables
#' @param sigma_t standard deviation for the treated group observations
#' @param sigma_c standard deviation for the control group observations
#'
#' @return a scalar value
#' @export
#'
#' @examples
a_fun <- function(delta, ht, hc, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c){
  
  if(missing(ht) | missing(hc)){
    
    hres <- bw_adpt(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
    ht <- hres$ht
    hc <- hres$hc
  }
  
  num_it1 <- K_fun(b = ht, C_pair = c(Cbar, Cj), X = Xt, mon_ind = mon_ind) / sigma_t^2
  num_it2 <- Cbar * Norm(Vplus(Xt, mon_ind)) - Cj * Norm(Vminus(Xt, mon_ind))
  num_it <- num_it1 * num_it2
  denom_it <- K_fun(b = ht, C_pair = c(Cbar, Cj), X = Xt, mon_ind = mon_ind) / sigma_t^2
  res_t <- sum(num_it) / sum(denom_it)
  
  num_ic1 <- K_fun(b = hc, C_pair = c(Cj, Cbar), X = Xc, mon_ind = mon_ind) / sigma_c^2
  num_ic2 <- Cj * Norm(Vplus(Xc, mon_ind)) - Cbar * Norm(Vminus(Xc, mon_ind))
  num_ic <- num_ic1 * num_ic2
  denom_ic <- K_fun(b = hc, gam_pair = c(1, 1), C_pair = c(Cj, Cbar), Xc, mon_ind) / sigma_c^2
  res_c <- sum(num_ic) / sum(denom_ic)
  
  res <- 0.5 * (res_t - res_c)
  return(res)
  
}

#' Worst-case bias of Lhat
#'
#' @param delta a nonegative scalar value: 
#' it can be left unspecified if `ht` and `hc` are specified
#' @param ht the modulus value for the treated observations: 
#' it can be left unspecified if `delta` is specified
#' @param hc the modulus value for the control observations: 
#' it can be left unspecified if `delta` is specified
#' @param Cj the smoothness parameter aiming to adapt to
#' @param Cbar the largest smoothness parameter
#' @param Xt data for the treated
#' @param Xc data for the control
#' @param mon_ind index for monotone variables
#' @param sigma_t standard deviation for the treated group observations
#' @param sigma_c standard deviation for the control group observations
#'
#' @return a scalar value
#' @export
#'
#' @examples
sup_bias_Lhat_RD <- function(delta, ht, hc, Cj, Cbar, Xt, Xc, mon_ind, 
                        sigma_t, sigma_c){
  
  if(missing(ht) | missing(hc)){
    
    hres <- bw_adpt(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
    ht <- hres$ht
    hc <- hres$hc
  }
  
  res1 <- 0.5 * (ht + hc)
  res2 <- a_fun(delta, ht, hc, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
  res3 <- -0.5 * (delta^2 / ht) / 
    sum(K_fun(b = ht, C_pair = c(Cbar, Cj), X = Xt, mon_ind = mon_ind) / sigma_t^2)
  
  res <- res1 + res2 + res3
  return(res)
}

#' standard deviation of Lhat
#'
#' @param delta a nonegative scalar value: 
#' it can be left unspecified if `ht` and `hc` are specified
#' @param ht the modulus value for the treated observations: 
#' it can be left unspecified if `delta` is specified
#' @param hc the modulus value for the control observations: 
#' it can be left unspecified if `delta` is specified
#' @param Cj the smoothness parameter aiming to adapt to
#' @param Cbar the largest smoothness parameter
#' @param Xt data for the treated
#' @param Xc data for the control
#' @param mon_ind index for monotone variables
#' @param sigma_t standard deviation for the treated group observations
#' @param sigma_c standard deviation for the control group observations
#'
#' @return
#' @export
#'
#' @examples
sd_Lhat_RD <- function(delta, ht, hc, Cj, Cbar, Xt, Xc, mon_ind, 
                      sigma_t, sigma_c){
 
  if(missing(ht) | missing(hc)){
    
    hres <- bw_adpt(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
    ht <- hres$ht
    hc <- hres$hc
  }
  
  res <- (delta / ht) / 
    sum(K_fun(b = ht, C_pair = c(Cbar, Cj), X = Xt, mon_ind = mon_ind) / sigma_t^2)
  
  return(res)
}

#' Lower CI for the RD parameter adapting to \eqn{C_j}
#'
#' @param delta a nonegative scalar value: 
#' it can be left unspecified if `ht` and `hc` are specified
#' @param ht the modulus value for the treated observations: 
#' it can be left unspecified if `delta` is specified
#' @param hc the modulus value for the control observations: 
#' it can be left unspecified if `delta` is specified
#' @param Cj the smoothness parameter aiming to adapt to
#' @param Cbar the largest smoothness parameter
#' @param Xt data for the treated
#' @param Xc data for the control
#' @param mon_ind index for monotone variables
#' @param sigma_t standard deviation for the treated group observations
#' @param sigma_c standard deviation for the control group observations
#' @param Yt outcome value for the treated group observations
#' @param Yc outcome value for the control group observations
#' @param tau 1 - (coverage probability)
#'
#' @return
#' @export
#'
#' @examples
c_hat_lower_RD <- function(delta, ht, hc, Cj, Cbar, Xt, Xc, mon_ind, 
                           sigma_t, sigma_c, Yt, Yc, tau) {
  
  lhat <- Lhat_fun_RD(delta, ht, hc, Cj, Cbar, Xt, Xc, mon_ind, 
                                  sigma_t, sigma_c, Yt, Yc)
  
  sup_bias <- sup_bias_Lhat_RD(delta, ht, hc, Cj, Cbar, Xt, Xc, mon_ind, 
                                           sigma_t, sigma_c)
  
  sd <- sd_Lhat_RD(delta, ht, hc, Cj, Cbar, Xt, Xc, mon_ind, 
                               sigma_t, sigma_c)
  
  res <- lhat - sup_bias - qnorm(1 - tau) * sd
  return(res)
}
