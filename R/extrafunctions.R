



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
  fits = huge::huge(t(B), method = method)
  solution = huge::huge.select(fits)
  return(Theta = as.matrix(solution$refit))
}

#' Calculate graph matching error
#'
#' This function counts the number of vertices that are incorrectly matched across two permutations, and divides
#' over the total number of vertices.
#'
#' @param Q1 Permutation matrix
#' @param Q2 Permutation matrix
#' @return Proportion of vertices with different matching.
#' @export
gm_error <- function(Q1, Q2) {
  n = ncol(Q1)
  sum(abs(Q1 - Q2)) / (2*n)
}
#' Calculate edge error
#'
#' This function counts the number edge disagreements between A and a permuted version of it, and
#' reports the proportion of edges that are different.
#' If the two oermutations are the same, then the result is zero.
#'
#' @param Q1 Permutation matrix
#' @param Q2 Permutation matrix
#' @param A adjacency matrix
#' @return Proportion of edges with different matching.
#' @export
edge_error <- function(Q1, Q2, A) {
  n = ncol(A)
  A = 1*(A!=0)
  edges = sum(A)
  sum(abs(Matrix::t(Q1) %*% A %*% (Q1) - Matrix::t(Q2) %*% A %*% (Q2))) / (2*edges)
}

#' Calculate errors in estimating edges of a graphical model
#'
#' This function reports the false positive rate, false negative rate, and proportion of errors for
#' estimating the edges of a graphical model.
#'
#' @param Theta A matrix where non-zeros indicate edges.
#' @param A Adjacency matrix of the graphical model
#' @param Q Alignment permutation matrix between A and Theta. The default is the identity matrix.
#' @return False positive rate (FPR), false negative rate (FNR) and total support error.
#' @export
graphical_model_errors <- function(Theta, A, Q = diag(ncol(A))) {
  n = ncol(A)
  A = 1*(A!=0)
  Theta = 1*(Theta!=0)
  diag(Theta) = 0
  Theta = Q %*% (1*(Theta != 0)) %*% Matrix::t(Q)
  diag(A) = 0
  FPR = sum( (Theta == 1) * (1-A)) / (sum((A==0)) - n)
  FNR = sum( (Theta==0) * A) / sum(A)
  support_Error = sum(Theta != A)/ (n*(n-1))
  return(c(FPR = FPR, FNR = FNR, support_error = support_Error))
}
