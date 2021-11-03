#only binomial implemented
penalized_pseudolikelihood <- function(B, family = "binomial", lambda, W, gamma = 1/100) {
  n = nrow(B)
  Theta_pseudo = Matrix(0, n, n)
  beta_pseudo = rep(0, n)
  for(k in 1:n) {
    fit_k <- glmnet::glmnet(x = t(B[-k,]), y = B[k,], family = family, #lambda = lambda,
                    penalty.factor = W[k,-k] + mean(W[k,-k])*gamma)
    # constant gamma added for numerical stability
    coefs =  glmnet::predict.glmnet(fit_k, s = lambda, type = "coefficients")
    Theta_pseudo[k, -k]= coefs[-1]
    beta_pseudo[k] = coefs[1]
  }
  return(list(beta = beta_pseudo, Theta = Theta_pseudo))
}

pseudolikelihood <- function(B, family = "binomial", U) {
  if(min(B) == -1 & length(unique(B)) == 2) # make binary to {-1, 1}
    B <- (B+1)/2
  n <- nrow(B)
  m <- ncol(B)
  Theta_pseudo <- Matrix(0, n, n)
  beta_pseudo <- rep(0, n)
  pseudo_vec <- c()
  for(k in 1:n) {
    dat <- data.frame(x = t(B[which(U[k,] !=0),,drop=FALSE]), y = B[k,])

    fit_k <- glm(formula = y~., data = dat, family = family)
    coefs <-  fit_k$coefficients

    Theta_pseudo[k, which(U[k,] != 0)] <- coefs[-1]
    beta_pseudo[k] <- coefs[1]
    if(family == "binomial") {
      pseudo_vec <- c(pseudo_vec, sum(log(1+exp(-predict(object = fit_k) * dat$y))))
    } else {
      pseudo_vec <- c(pseudo_vec, -(dat$y - predict(object = fit_k))^2)
    }
  }
  return(list(beta = beta_pseudo, Theta = Theta_pseudo, pseudolik = sum(pseudo_vec)/m))
}
