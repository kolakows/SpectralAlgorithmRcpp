Laplacian_eigen <- function(G,k){
  stopifnot(k>1)
  D <- diag(rowSums(G))
  L <- D-G
  V <- eigen(L,symmetric = TRUE)[[2]]
  n <- ncol(G)
  V[,(n-k):(n-1)]
}

spectral_clustering <- function(X, k, M){
  library(Rcpp)
  #load Mnn and Mnn_graph functions from .cpp file
  sourceCpp("./spectral_aux.cpp",rebuild=TRUE)
  S <- Mnn(as.matrix(X), M)
  G <- Mnn_graph(S)
  E <- Laplacian_eigen(G,k)
  kmeans(E,k,iter.max = 30,nstart=5)
}

