#' Bipartite to unipartite matching via collapsing the bipartite network
#'
#' This function performs graph matching between a bipartite and a unipartite graphs that are assumed to
#' share a common set of vertices. The method performs an alignment of the edges of the unipartite graphs
#' with the edges of a collapsed bipartite graph into a unipartite, following the collapsed representation
#' selected.
#'
#' @param A Adjacency matrix of the unipartite graph.
#' @param B Bipartite graph incidence matrix, where rows correspond to the vertices in common
#' with the unipartite graph.
#' @param collapsed_method Collapsed representation for the bipartite graph. The options are one-mode projection (\texttt{omp}),
#' correlation (\texttt{corr}), covariance (\texttt{cov}), graphical lasso (\texttt{glasso}) or nodewise regression (\texttt{mb}).
#' @param seeds If some vertices have known correspondence, a vector containing the indexes of these
#' vertices can be passed through this parameter, and the corresponding rows of A and B are assumed to be aligned.
#' The algorithm will then match the remaining vertices. The default is NULL if no seeds are available.
#' @return A permutation Q that approximately minimizes the difference between the edges of A and the collapsed B.
#' @export
bipartite_matching_collapsed <- function(A, B, collapsed_method = c("omp", "corr", "cov",
                                                                              "glasso", "mb"),
                                         seeds = NULL){

  collapsed_method <- match.arg(collapsed_method, c("omp", "corr", "cov",
                                                    "glasso", "mb"))
  m <- ncol(B)
  n <- nrow(B)
  Bcollapsed <- switch (collapsed_method,
    omp = tcrossprod(B)/m,
    corr = cor(Matrix::t(B)),
    cov = cov(Matrix::t(B)),
    glasso = graphical_model_estimation(B, method = "glasso"),
    mb = graphical_model_estimation(B, method = "mb"))

  if(sum(abs(Bcollapsed)) == 0) {
    return(list(Q = diag(ncol(A))[sample(1:ncol(A), ncol(A)),]))
  }


  GM = iGraphMatch::graph_match_FW(A = A, B = Bcollapsed,
                                   seeds = seeds,
                                   start = "bari")
  Q <- GM$P
  return(list(Q = Q, Bcollapsed = Bcollapsed))
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
  Theta_hat = suppressWarnings(glasso::glasso(S, rho=0, zero = which(Q_1==1, arr.ind = T)))
  f_t = log(det(Theta_hat$wi)) - sum(S*Theta_hat$wi)
  return(list(Q = Q, Theta = Theta_hat$wi, f = f_t))
}
