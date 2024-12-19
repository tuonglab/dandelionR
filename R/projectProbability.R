#' Project Probabilities from Markov Chain to Pseudobulks
#'
#' This function projects probabilities calculated from a Markov chain onto each pseudobulk based on a diffusion distance matrix.
#'
#' @param diffusionmap diffusion map, used to reconstruct diffustion distance matrix
#' @param waypoints Integer vector. Indices of the waypoints used in the Markov chain.
#' @param probabilities Numeric vector. Probabilities associated with the waypoints, calculated from the Markov chain.
#' @param t Numeric. The diffusion time to be used in the projection.
#' @importFrom destiny eigenvectors eigenvalues
#' @importFrom stats sd
#' @return each pseudobulk's probabilites
projectProbability <- function(diffusionmap, waypoints, probabilities, t = 1) {
    message("Project probabilites from waypoints to each pseudobulk...")
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
    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            diff <- eigenvectors[i, seq_len(K)] - eigenvectors[j, seq_len(K)]
            D_diffusion[i, j] <- sqrt(sum(lambda_t * (diff^2)))
        }
    }
    # `D_diffusion` is now the diffusion distance matrix
    Dif <- D_diffusion[, waypoints]
    D_ravel <- as.vector(Dif)
    requireNamespace("stats")
    sdv <- stats::sd(D_ravel) * 1.06 * (length(D_ravel)^(-1 / 5))
    W <- exp(-0.5 * ((Dif / sdv)^2))
    W <- W / apply(W, 1, sum)
    prob <- W %*% probabilities
    return(prob)
}
