#' @export
bipartite_matching_collapsed <- function(A, B_collapsed, S=NULL, seeds = NULL,
                                         similarity = FALSE) {
  if(sum(abs(B_collapsed)) == 0) {
    return(list(Q = diag(ncol(A))[sample(1:ncol(A), ncol(A)),]))
  }

  if(similarity) {
    embedA <- ase(A)
    embedB <- ase(B_collapsed)
    d = max(ncol(embedA), ncol(embedB))
    embedA <- ase(A, d)
    embedB <- ase(B_collapsed, d)
    simMatrix <- tcrossprod(embedA, embedB)
  } else {
    simMatrix <- NULL
  }

  n = ncol(A)
  m = ncol(B_collapsed)

  GM = iGraphMatch::graph_match_FW(A = A, B = B_collapsed,
                                   seeds = seeds,
                                   similarity = simMatrix,
                                   start = "bari")
  Q = GM$P
  Q_1 = as.matrix(1 - t(Q) %*% A %*% Q)
  diag(Q_1) = 0
  if(is.null(S)){
    Theta_hat = NULL
    f_t = NULL
  } else{
    require(glasso)
    Theta_hat = suppressWarnings(glasso(S, rho=0, zero = which(Q_1==1, arr.ind = T)))
    f_t = log(det(Theta_hat$wi)) - sum(S*Theta_hat$wi)
  }
  return(list(Q = Q, Theta = Theta_hat$wi, f = f_t))
}




bipartite_matching_collapsed_starts <- function(A, B_collapsed, S, seeds = NULL, nstarts = 10) {
  n = ncol(A)
  m = ncol(B_collapsed)
  GMlist = list()
  for(i in 1:nstarts)
    GMlist[[i]] = iGraphMatch::graph_match_FW(A = A, B = B_collapsed, seeds = seeds,
                                              start = "bari", max_iter  = 100)
  Asum = sum(abs(A))
  B_col = B_collapsed / sum(abs(B_collapsed))*Asum
  loss_function <- sapply(GMlist, function(GM) {
    P = GM$P
    sum(abs(A - P %*% B_col %*% t(P))^2)
  })
  GM = GMlist[[which.min(loss_function)]]
  Q = GM$P
  Q_1 = as.matrix(1 - t(Q) %*% A %*% Q)
  diag(Q_1) = 0
  Theta_hat = suppressWarnings(glasso(S, rho=0, zero = which(Q_1==1, arr.ind = T)))
  f_t = log(det(Theta_hat$wi)) - sum(S*Theta_hat$wi)
  return(list(Q = Q, Theta = Theta_hat$wi, f = f_t))
}
