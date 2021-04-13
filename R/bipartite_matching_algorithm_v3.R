# Input: A, B, grid of lambda values
# Ouput: Inverse covariance, Alignment matrix

#' Bipartite to unipartite matching via penalized inverse covariance estimation
#'
#' This function performs graph matching between a bipartite and a unipartite graph that are assumed to
#' share a common set of vertices. The method performs an alignment of the non-zero entries of the inverse covariance matrix
#' of the bipartite graph, and the nonzero entries of the unipartite adjacency matrix.
#'
#' @param A Adjacency matrix of the unipartite graph.
#' @param B Bipartite graph incidence matrix, where rows correspond to the vertices in common
#' with the unipartite graph.
#' @param Q_true A permutation matrix with the true solution, if known. When the true solution is unknown, the default is NULL.
#' @param lambdas Penalty values for the optimization algorithm.
#' @param MAX_ITER Maximum number of iterations.
#' @param verbose If TRUE, the method prints an output after each iteration.
#' @param covariance Wheter to use covariance or correlation matrix.
#' @param seeds If some vertices have known correspondence, a vector containing the indexes of these
#' vertices can be passed through this parameter, and the corresponding rows of A and B are assumed to be aligned.
#' The algorithm will then match the remaining vertices. The default is NULL if no seeds are available.
#' @return
#' @export
bipartite_matching_icov <- function(A, B,
                               Q_true = NULL,
                               lambdas = 10^(seq(-2, -1, length.out = 5)),
                               MAX_ITER = 30,
                               verbose = FALSE,
                               covariance = TRUE,
                               seeds = NULL) {
  # algorithm parameters
  TOL <- 1e-4
  l_factor <- 1

  n <- nrow(B)
  m <- ncol(B)

  if(covariance) {
    S <- cov(Matrix::t(B))
  }else {
    S <- cor(Matrix::t(B))
  }

  # check if zero rows in B
  zerorowsB <- sum(rowSums(B!=0)==0) > 0
  if(zerorowsB) whichzerorowsB <- which(rowSums(B!=0)==0)

  # penalty values
  Q_list_lambda <- list()
  f_list_lambda <- list()

  # report baseline loglikelihood if true solution is known
  if(!is.null(Q_true)) {
    Q_1 <- as.matrix(1 - Matrix::t(Q_true) %*% A %*% Q_true)
    diag(Q_1) <- 0
    Theta_hat <- suppressWarnings(glasso::glasso(S, rho=0, zero = which(Q_1==1, arr.ind = T)))
    if(is.nan(max(Theta_hat$wi))) {
      Theta_hat <- suppressWarnings(glasso::glasso(S, rho=1e-6, zero = which(Q_1==1, arr.ind = T)))
      if(verbose)
        cat(sprintf("---Baseline penalized loglikelihood (rho=1e-6): %.4f\n",
            log(det(Theta_hat$wi)) - sum(S*Theta_hat$wi)))
    } else {
      if(verbose)
        cat(sprintf("---Baseline loglikelihood: %.4f\n", log(det(Theta_hat$wi)) - sum(S*Theta_hat$wi)))
    }
  }

  i = 1
  for(lambda in lambdas) {
    if(verbose ==TRUE) cat("-Running lambda=", lambda, "\n")

    # run algorithm --------------------------------------------
    # initialize
    iter <- 0
    crit <- Inf
    D <- matrix(1/n, n, n)
    f_0 <- -Inf
    Q_list <- list()
    f_list <- list()

    while(iter <= MAX_ITER & crit > TOL) {
      iter <- iter + 1
      # make penalty matrix
      P <- Matrix::t(D) %*% A %*% (D)
      diag(P) <- 1
      P <- 1- P
      P <- as.matrix((P + Matrix::t(P))/2) # make sure it is symmetric

      # run graphical lasso =================================
      if(zerorowsB) { # delete zero rows of B if they exist
        GLsub <- glasso::glasso(S[-whichzerorowsB, -whichzerorowsB],
                       rho = lambda*(P[-whichzerorowsB, -whichzerorowsB]))
        # penalize all entries if solution is not found
        if(is.nan(sum(GLsub$wi[-whichzerorowsB, -whichzerorowsB])))
          GLsub <- glasso::glasso(S[-whichzerorowsB, -whichzerorowsB],
                         rho = lambda*(P[-whichzerorowsB, -whichzerorowsB]) + lambda/10)
        # absolute inverse covariance
        U <- matrix(0, n, n)
        U[-whichzerorowsB, -whichzerorowsB] <- abs(GLsub$wi)
      } else{ # B has no zero rows
        GL <- glasso::glasso(S, rho = lambda*(P))
        if(is.nan(sum(GL$wi)))
          GL <- glasso::glasso(S, rho = lambda*(P) + lambda/10)
        # inverse covariance
        U <- abs(GL$wi)
      }
      diag(U) <- 0
      if(zerorowsB) {
        U[whichzerorowsB,]  <-  U[,whichzerorowsB] = 0
      }
      sparse  <-  sum(U!=0)/(n^2-n) # calculate %off-diagonal nonzeros

      # run graph matching ====================================
      GM <- iGraphMatch::graph_match_FW(A = A, B = U,
                                       start = "bari", seeds = seeds)
      # graph matching solution
      D <- as.matrix(GM$D)
      # compute objective function value
      Q <- GM$P
      Q_1 <- as.matrix(1 - Matrix::t(Q) %*% A %*% (Q))
      diag(Q_1) <- 0
      Theta_hat <- suppressWarnings(glasso::glasso(S, rho=0, zero = which(Q_1==1, arr.ind = T)))
      if(is.nan(max(Theta_hat$wi)))
        Theta_hat <- suppressWarnings(glasso::glasso(S, rho=1e-6, zero = which(Q_1==1, arr.ind = T)))
      if(zerorowsB){
        f_t <- log(det(Theta_hat$wi[-whichzerorowsB, -whichzerorowsB])) -
          sum((S*Theta_hat$wi)[-whichzerorowsB,-whichzerorowsB])
      } else{
        f_t <- log(det(Theta_hat$wi)) - sum(S*Theta_hat$wi)
      }
      lambda <- lambda * l_factor
      crit <- abs(f_0 -f_t)
      f_0 <- f_t
      # save parameters
      Q_list[[iter]] <- Q
      f_list[[iter]] <- f_t

      # Print iteration results (% nonzeros, % incorrect matched vertices)
      if(verbose == TRUE) {
        if(is.null(Q_true)) {
          cat(sprintf(c("--- %d\tloglik_t=%.4f\t%%nonzeros=%.4f\n", iter, f_t, sparse)))
        } else{
          if(is.null(seeds)) {
            cat(sprintf("---%d\tloglik_t= %0.4f\t%%nonzeros= %0.4f\tMatch.%%error= %0.4f\n",
                        iter, f_t, sparse,100*sum(abs(Q_true - D))/(2*(n)), "\n" ))
          } else{
            cat(sprintf("---%d\tloglik_t= %0.4f\t%%nonzeros= %0.4f\tMatch.%%error= %0.4f\n",
                        iter, f_t, sparse,
                        100*sum(abs(Q_true[-seeds, -seeds] - D[-seeds, -seeds]))/(2*(n-length(seeds))), "\n" ))
          }

        }
      }
    }
    # choose best solution
    index <- which.max(unlist(f_list))
    Q_list_lambda[[i]] <- Q_list[[index]]
    f_list_lambda[[i]] <- f_list[[index]]
    i <- i + 1
  }

  # overall best across lambdas
  index <- which.max(unlist(f_list_lambda))
  Q_best <- Q_list_lambda[[index]]
  f_best <- f_list_lambda[[index]]
  Q_1 <- as.matrix(1 - Matrix::t(Q_best) %*% A %*% Q_best)
  diag(Q_1) <- 0
  Theta_hat <- suppressWarnings(glasso::glasso(S, rho=0, zero = which(Q_1==1, arr.ind = T)))

  return(list(Q = Q_best, Theta = Theta_hat$wi, f = f_best))
}


bipartite_matching_fast <- function(A, B, Q_true = NULL,
                                    method = "glasso",
                                    verbose = FALSE,
                                    covariance = TRUE,
                                    seeds = NULL) {
  n = nrow(B)
  m = ncol(B)
  if(covariance) {
    S = cov(Matrix::t(B))
  }else {
    S = cor(Matrix::t(B))
  }
  hg = huge(Matrix::t(B), method = method)
  # penalty values
  Q_list_lambda = list()
  f_list_lambda = list()
  i = 1
  # baseline
  if(!is.null(Q_true)) {
    Q_1 = as.matrix(1 -Matrix::t(Q_true) %*% A %*% Q_true)
    diag(Q_1) = 0
    Theta_hat = suppressWarnings(glasso::glasso(S, rho=0, zero = which(Q_1==1, arr.ind = T)))
    if(is.nan(max(Theta_hat$wi)))
      Theta_hat = suppressWarnings(glasso::glasso(S, rho=1e-4, zero = which(Q_1==1, arr.ind = T)))
    print(sprintf("---Baseline logdet: %.4f", log(det(Theta_hat$wi)) - sum(S*Theta_hat$wi)))
  }
  # screen
  nonzeros = sapply(hg$path, sum)
  if(sum(nonzeros) == 0) stop("empty graphs were estimated")
  for(i in which(nonzeros>0)) {
    if(verbose ==TRUE)
      cat("-Running lambda ", i, "\n")
    GM = iGraphMatch::graph_match_FW(A = A, B = hg$path[[i]],
                                     start = "bari", seeds = seeds)
    # matching solution
    D = as.matrix(GM$D)

    # objective function to calculate convergence and fit
    Q = GM$P
    Q_1 = as.matrix(1 - Matrix::t(Q) %*% A %*% Q)
    diag(Q_1) = 0
    Theta_hat = suppressWarnings(glasso::glasso(S, rho=0, zero = which(Q_1==1, arr.ind = T)))
    if(is.nan(max(Theta_hat$wi)))
      Theta_hat = suppressWarnings(glasso::glasso(S, rho=1e-3, zero = which(Q_1==1, arr.ind = T)))
    f_t = log(det(Theta_hat$wi)) - sum(S*Theta_hat$wi)

    # choose best solution
    Q_list_lambda[[i]] = Q
    f_list_lambda[[i]] = f_t
  }
  # overall best
  index = which.max(unlist(f_list_lambda))
  Q_best = Q_list_lambda[[index]]
  f_best = f_list_lambda[[index]]
  Q_1 = as.matrix(1 - Matrix::t(Q_best) %*% A %*% Q_best)
  diag(Q_1) = 0
  Theta_hat = suppressWarnings(glasso::glasso(S, rho=0, zero = which(Q_1==1, arr.ind = T)))

  return(list(Q = Q_best, Theta = Theta_hat$wi, f = f_best))
}





