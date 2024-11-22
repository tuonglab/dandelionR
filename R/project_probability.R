#' project_probabilty
#'
#' project the probabilities from Markov chain to each pseudobulk
#' @param diffusionmap diffusion map, used to reconstruct diffustion distance matrix
#' @param waypoints index of waypoints
#' @param probabilities waypoints' probabilities, result from markov chain
#' @param t diffusion time
#' @return each pseudobulk's probabilites
project_probability <- function(diffusionmap, waypoints, probabilities, t = 1) {
  message("Project probabilites from waypoints to each pseudobulk...")
  # Extract eigenvalues and eigenvectors from the DiffusionMap
  requireNamespace("destiny")
  eigenvalues <- destiny::eigenvalues(diffusionmap)  # vector of eigenvalues
  eigenvectors <- destiny::eigenvectors(diffusionmap)  # matrix of eigenvectors (each column is an eigenvector)
  # Set diffusion time `t` and the number of components `K` to use
  t <- 1  # diffusion time
  K <- min(length(eigenvalues), ncol(eigenvectors))  # use available components
  # Precompute squared eigenvalues scaled by diffusion time
  lambda_t <- (eigenvalues[1:K])^(2 * t)
  # Initialize an empty matrix for the diffusion distances
  n <- nrow(eigenvectors)
  D_diffusion <- matrix(0, nrow = n, ncol = n)
  # Calculate the pairwise diffusion distance
  for (i in 1:n) {
    for (j in 1:n) {
      diff <- eigenvectors[i, 1:K] - eigenvectors[j, 1:K]
      D_diffusion[i, j] <- sqrt(sum(lambda_t * (diff^2)))
    }
  }
  # `D_diffusion` is now the diffusion distance matrix
  Dif <- D_diffusion[, waypoints]
  D_ravel <- as.vector(Dif)
  requireNamespace("stats")
  sdv <- stats::sd(D_ravel) * 1.06 * (length(D_ravel)^(-1/5))
  W = exp(-0.5 * ((Dif/sdv)^2))
  W = W/apply(W, 1, sum)
  prob <- W %*% probabilities
  return(prob)
}