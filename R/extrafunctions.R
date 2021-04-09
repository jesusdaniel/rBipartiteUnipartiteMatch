



loss_glasso_theta <- function(Q, A, B) {
  require(glasso)
  Q_1 = as.matrix(1 - t(Q) %*% A %*% (Q))
  diag(Q_1) = 0
  S = cov(t(B))
  zerorowsB = sum(rowSums(B!=0)==0)>0
  if(zerorowsB) whichzerorowsB = which(rowSums(B!=0)==0)


  if(zerorowsB){
    Theta_hat = suppressWarnings(glasso(S[-whichzerorowsB, -whichzerorowsB],
                                        rho=1e-6,
                                        zero = which(Q_1[-whichzerorowsB, -whichzerorowsB]==1, arr.ind = T)))
    f_t = log(det(Theta_hat$wi)) -
      sum((S[-whichzerorowsB, -whichzerorowsB]*Theta_hat$wi))
  } else{
    Theta_hat = suppressWarnings(glasso(S,
                                        rho=1e-4, zero = which(Q_1==1, arr.ind = T)))
    f_t = log(det(Theta_hat$wi)) - sum(S*Theta_hat$wi)
  }
  return(f_t)
}


# methods = c("mb", "glasso")
graphical_model_estimation <- function(B, method = "mb") {
  fits = huge(t(B), method = method)
  solution = huge.select(fits)
  return(Theta = as.matrix(solution$refit))
}
