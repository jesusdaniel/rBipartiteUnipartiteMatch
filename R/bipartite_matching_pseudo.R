
# Input: A, B, grid of lambda values
# Ouput: Inverse covariance, Alignment matrix
#' Bipartite to unipartite matching via penalized pseudolikelihood.
#'
#' This function performs graph matching between a bipartite and a unipartite graphs that are assumed to
#' share a common set of vertices. The method performs an alignment of the non-zero parameters of a 
#' graphical model for the bipartite graph and the nonzero entries of the unipartite adjacency matrix
#' via pseudolikelihood maximization.
#'
#' @param A Adjacency matrix of the unipartite graph.
#' @param B Bipartite graph incidence matrix, where rows correspond to the vertices in common
#' with the unipartite graph.
#' @param Q_true A permutation matrix with the true solution, if known. When the true solution is unknown, the default is NULL.
#' @param lambdas Penalty values for the optimization algorithm.
#' @param MAX_ITER Maximum number of iterations.
#' @param verbose If TRUE, the method prints an output after each iteration.
#' @param seeds If some vertices have known correspondence, a vector containing the indexes of these
#' vertices can be passed through this parameter, and the corresponding rows of A and B are assumed to be aligned.
#' The algorithm will then match the remaining vertices. The default is NULL if no seeds are available.
#' @param gamma Additional penalty for all variables. A positive small constant is needed when the sample size
#' @param similarity Whether to use a similarity matrix in graph matching (experimental parameter).
#' @return
#' @export
bipartite_matching_pseudolikelihood <- function(A, B, Q_true = NULL,
                                                family = "binomial",
                                                lambdas = 10^(-2:0),
                                                MAX_ITER = 20, verbose = FALSE,
                                                seeds = NULL,
                                                gamma = 0,
                                                similarity = FALSE) {
  # algorithm parameters
  TOL <- 1e-4
  l_factor <- 1
  
  n <- nrow(B)
  m <- ncol(B)
  
  # penalty values
  Q_list_lambda <- list()
  f_list_lambda <- list()
  
  
  if(similarity) {
    embedA <- ase(A)
    dimembed <- ncol(embedA)
    simMatrix <- matrix(0, n, n)
  } else {
    simMatrix <- NULL
  }
  
  # check if zero rows in B
  zerorowsB <- sum(rowSums(B!=0)==0) > 0
  if(zerorowsB) whichzerorowsB <- which(rowSums(B!=0)==0)
  
  # report baseline loglikelihood if true solution is known
  if(!is.null(Q_true)) {
    Q_1 <- as.matrix(t(Q_true) %*% A %*% Q_true)
    diag(Q_1) <- 0
    ground_Truth <- pseudolikelihood(B, U = Q_1, family = family)
    if(is.null(ground_Truth))
      error("Error: too many non-zeros in A, not enough samples.")
    if(verbose ==TRUE) cat("---Baseline pseudolikelihood:", ground_Truth$pseudolik, "\n")
  }
  
  i <- 1
  for(lambda in lambdas) {
    if(verbose ==TRUE) cat("-Running lambda=", lambda, "\n")
    
    # run algorithm ---------------------------------------
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
      P <- t(D) %*% A %*% D
      diag(P) <- 1
      P <- 1- P
      P <- as.matrix((P + t(P))/2) # make sure it is symmetric
      
      # run penalized pseudolikelihood optimization =================================
      PL <- penalized_pseudolikelihood(B, lambda = lambda, W = P, family = family, 
                                      gamma = gamma)
      if(is.null(PL))
        error("Error: too many non-zeros in A, not enough samples. Increase gamma parameter.")
      # inverse covariance
      U <- (abs(PL$Theta))
      diag(U) <- 0
      
      sparse <- sum(U!=0)/(n^2-n)
      
      # run graph matching ========================================================
      if(similarity) {
        embedB <- ase(U + t(U), dimembed) / sqrt(2)
        simMatrix <- tcrossprod(embedA, embedB)
      }
      GM = iGraphMatch::graph_match_FW(A = A, B = (U + t(U))/2,
                                       similarity = simMatrix,
                                       start = "bari", seeds = seeds)
      # matching solution
      D <- as.matrix(GM$D)
      
      # objective function to calculate convergence and fit
      Q <- GM$P
      Q_1 <- as.matrix(t(Q) %*% A %*% Q)
      diag(Q_1) <- 0
      Pseudo_unpen <- pseudolikelihood(B, family = family, U = Q_1)
      if(is.null(Pseudo_unpen))
        error("Error: too many non-zeros in A, not enough samples.")

      f_t <- Pseudo_unpen$pseudolik
      lambda <- lambda * l_factor
      crit <- abs(f_0 -f_t)
      f_0 <- f_t
      # save parameters
      Q_list[[iter]] <- Q
      f_list[[iter]] <- f_t
      if(verbose == TRUE)
        cat(c("---", iter, "\t", "pseudologlik_t=", f_t, "\t %nonzeros=", sparse,
              ifelse(is.null(Q_true), "", c("Match.%error==")),
              ifelse(is.null(Q_true), "",ifelse(is.null(seeds), c(sum(abs(Q_true - D))/(2*(n))),
                                                c(sum(abs(Q_true[-seeds, -seeds] - D[-seeds, -seeds]))/(2*(n-length(seeds)))))),
              "\n"))
    }
    # choose best solution
    index <- which.max(unlist(f_list))
    Q_list_lambda[[i]] <- Q_list[[index]]
    f_list_lambda[[i]] <- f_list[[index]]
    i <- i + 1
  }
  # overall best
  index <- which.max(unlist(f_list_lambda))
  Q_best <- Q_list_lambda[[index]]
  f_best <- f_list_lambda[[index]]
  Q_1 <- as.matrix(Q_best %*% A %*% t(Q_best))
  diag(Q_1) <- 0
  Pseudo_unpen <- pseudolikelihood(B, family = family, U = Q_1)
  Theta_hat <- Pseudo_unpen$Theta
  
  return(list(Q = Q_best, Theta_hat=Theta_hat,f = f_best))
}


