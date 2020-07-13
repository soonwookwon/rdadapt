pos <- function(X) {
  return(pmax(X, 0))
}

Vplus <- function(X, mon_ind) {
  
  Xplus <- X
  Xplus[, mon_ind] <- pos(X[, mon_ind])
  
  return(Xplus)
}

Vminus <- function(X, mon_ind) {
  
  return(Vplus(-X, mon_ind))
}

K_fun <- function(b, gam_pair = c(1, 1), C_pair, X, mon_ind, swap = FALSE){
# gam_pair is for adaptation to Holder cofficients
  
  if (swap) {
    gam_pair <- gam_pair[2:1]
    C_pair <- C_pair[2:1]
  }
  
  K <- pos(1 - (C_pair[1] / b) * Norm(Vplus(X, mon_ind))^gam_pair[1] -
             (C_pair[2] / b) * Norm(Vminus(X, mon_ind))^gam_pair[2])
  
  return(K)
}
