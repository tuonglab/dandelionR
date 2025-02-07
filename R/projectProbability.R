#' Project Probabilities from Markov Chain to Pseudobulks
#'
#' This function projects probabilities calculated from a Markov chain onto each pseudobulk based on a diffusion distance matrix.
#'
#' @param diffusionmap diffusion map, used to reconstruct diffustion distance matrix
#' @param waypoints Integer vector. Indices of the waypoints used in the Markov chain.
#' @param probabilities Numeric vector. Probabilities associated with the waypoints, calculated from the Markov chain.
#' @param t Numeric. The diffusion time to be used in the projection.
#' @param verbose Boolean, whether to print messages/warnings.
#' @importFrom destiny eigenvectors eigenvalues
#' @importFrom stats sd
#' @return each pseudobulk's probabilites
projectProbability <- function(diffusionmap, waypoints, probabilities, t = 1, verbose = TRUE) {
    if (verbose) {
        message("Project probabilites from waypoints to each pseudobulk...")
    }
    # Extract eigenvalues and eigenvectors from the DiffusionMap
    eigenvalues <- eigenvalues(diffusionmap) # vector of eigenvalues
    eigenvectors <- eigenvectors(diffusionmap) # matrix of eigenvectors (each column is an eigenvector)
    # Set diffusion time `t` and the number of components `K` to use
    t <- 1 # diffusion time
    K <- min(length(eigenvalues), ncol(eigenvectors)) # use available components
    # Precompute squared eigenvalues scaled by diffusion time
    lambda_t <- (eigenvalues[seq_len(K)])^(2 * t)
    # Initialize an empty matrix for the diffusion distances
    n <- nrow(eigenvectors)
    D_diffusion <- matrix(0, nrow = n, ncol = n)
    # Calculate the pairwise diffusion distance
    D_diffusion <- Reduce(function(dfm_j,j) {
      dfm_j <- Reduce(function(dfm_i, i){
        calDif(dfm, i, j = j, eigenvectors = eigenvectors, lambda_t = lambda_t, K = K)
      }, seq_len(n), init = dfm_j)
      dfm_j
    }, seq_len(n), init = D_diffusion)
    # `D_diffusion` is now the diffusion distance matrix
    Dif <- D_diffusion[, waypoints]
    D_ravel <- as.vector(Dif)
    sdv <- sd(D_ravel) * 1.06 * (length(D_ravel)^(-1 / 5))
    W <- exp(-0.5 * ((Dif / sdv)^2))
    W <- W / apply(W, 1, sum)
    prob <- W %*% probabilities
    return(prob)
}


#' function help to reconstruct diffustion distance using Reduce
#' 
#' @param dfm an nrow(eigenvectors) x nrow(eigenvectors) matrix
#'  need to be filled with diffusion distance with iteration
#' @param i integer the row selected in this iteration
#' @param j integer the col selected in this iteration
#' @param eigenvector numeric vector, the eigenvector from diffusion map
#' @param lambda_t eigenvalues to the power of t(diffusion time)
#' @param K The number of the eigenvectors to be used in calculation
#' @keywords internal
#' @return updated diffusion distance matrix after one iteration
.calDif <- function(dfm, i, j, eigenvectors, lambda_t, K)
{
  diff <- eigenvectors[i, seq_len(K)] - eigenvectors[j, seq_len(K)]
  D_diffusion[i, j] <- sqrt(sum(lambda_t * (diff^2)))
  D_diffusion
}
