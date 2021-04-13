##############################################
# GRAPH MATCHING BETWEEN BIPARTITE AND UNIPARTITE NETWORKS
#
# This script generates a pair of random networks, and applies
# different matching algorithms to align the vertices
# The unipartite network is generated using a Erdos-Renyi model, and
# the bipartite network follows the Ising model.

library(rBipartiteUnipartiteMatch)

set.seed(111)

#####################################################
# 1) Sample pair of graphs
# Number of vertices in the first class (common to both graphs)
n <- 100
# Number of vertices in the second class (bipartite graph)
m <- 1000
# Edge probability in unipartite graph A
p <- 0.05
# Sample unipartite graph
G_A <- igraph::sample_gnp(n, p)  
A <- igraph::get.adjacency(G_A)

# Sample bipartite graph
Q <- diag(n) # Correct permutation solution
W <- t(Q) %*% A %*% Q
mus <- rep(0, n)
max_weight <- 0.2
B <- t(IsingSampler::IsingSampler(n = m, graph = as.matrix(W), thresholds = mus, 
                     responses = c(-1,1),
                     beta = max_weight, method = "MH"))
B <- (B+1)/2

#####################################################
# 2) Bipartite graph matching
# Algorithm parameters
lambdas_inv <- 10^(seq(-2.5, -1, length.out = 10))
lambdas_pseudo <- 10^(seq(-2, -0.5, length.out = 10))
MAX_ITER_inv <- 20
MAX_ITER_pseudo <- 20

# Methods
BipGM1 <- bipartite_matching_icov(A, B, verbose = T, lambdas = lambdas_inv, 
                                 MAX_ITER = MAX_ITER_inv, seeds = NULL, Q_true = Q)
Bip_pseudoGM <- bipartite_matching_pseudolikelihood(A, B, verbose = T, 
                                                   lambdas = lambdas_pseudo, 
                                                   MAX_ITER = MAX_ITER_pseudo, 
                                                   Q_true = Q)

# Evaluate results
bipartite_gm_errors <- c(gm_error(BipGM1$Q, Q), gm_error(Bip_pseudoGM$Q, Q))
bipartite_edge_errors <- c(edge_error(BipGM1$Q, Q, A), edge_error(Bip_pseudoGM$Q, Q, A))
bipartite_graphmodel_errors <- cbind(graphical_model_errors(BipGM1$Q, A, Q), graphical_model_errors(Bip_pseudoGM$Q, A, Q))

#####################################################
# 3) Bipartite graph matching using collapsed bipartite networks
collapsed_results <- lapply(c("omp", "cov", "glasso", "mb"), function(method)
  bipartite_matching_collapsed(A = A, B, collapsed_method = method))

collapsed_gm_errors <- sapply(collapsed_results, function(method)
  gm_error(method$Q, Q))
collapsed_edge_errors <- sapply(collapsed_results, function(method)
  edge_error(method$Q, Q, A))
collapsed_graphicalmodel_errors <- sapply(collapsed_results, function(method)
  graphical_model_errors(method$Q, A, Q))

# Compare results
result_gm_errors <- c(bipartite_gm_errors, collapsed_gm_errors)
result_edge_errors <- c(bipartite_edge_errors, collapsed_edge_errors)
result_graphicalmodel_errors <- cbind(bipartite_graphmodel_errors, collapsed_graphicalmodel_errors)

results <- rbind(result_gm_errors,result_edge_errors,result_graphicalmodel_errors )

colnames(results) <- c("B-InvCov", "B-Pseudo", "C-OMP", "C-Cov", "C-GLasso", "C-M&B")

results