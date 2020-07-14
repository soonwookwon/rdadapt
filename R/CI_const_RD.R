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
  
  num_it <- K_fun(b = ht, gam_pair = c(1, 1), C_pair = c(Cbar, Cj), Xt, mon_ind) * Yt / sigma_t
  denom_it <- K_fun(b = ht, gam_pair = c(1, 1), C_pair = c(Cbar, Cj), Xt, mon_ind) / sigma_t
  res_t <- sum(num_it) / sum(denom_it)
  
  num_ic <- K_fun(b = hc, gam_pair = c(1, 1), C_pair = c(Cj, Cbar), Xc, mon_ind) * Yc / sigma_c
  denom_ic <- K_fun(b = hc, gam_pair = c(1, 1), C_pair = c(Cj, Cbar), Xc, mon_ind) / sigma_c
  res_c <- sum(num_ic) / sum(denom_ic)
  
  return(res_t - res_c)
}