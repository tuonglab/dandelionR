#' Project Probabilities from Markov Chain to Pseudobulks
#'
#' This function projects probabilities calculated from a Markov chain onto
#'  each pseudobulk based on a diffusion distance matrix.
#'
#' @param diffusionmap diffusion map, used to reconstruct
#'  diffustion distance matrix
#' @param waypoints Integer vector. Indices of the waypoints used in
#'  the Markov chain.
#' @param probabilities Numeric vector. Probabilities associated with
#'  the waypoints, calculated from the Markov chain.
#' @param t Numeric. The diffusion time to be used in the projection.
#' @param verbose Boolean, whether to print messages/warnings.
#' @importFrom destiny eigenvectors eigenvalues
#' @importFrom stats sd
#' @return each pseudobulk's probabilites
projectProbability <- function(
    diffusionmap, waypoints, probabilities, t = 1,
    verbose = TRUE) {
    if (verbose) {
        message("Project probabilites from waypoints to each pseudobulk...")
    }
    # Extract eigenvalues and eigenvectors from the DiffusionMap
    eigenvalues <- eigenvalues(diffusionmap) # vector of eigenvalues
    # matrix of eigenvectors (each column is an eigenvector)
    eigenvectors <- eigenvectors(diffusionmap)
    # Set diffusion time `t` and the number of components `K` to use
    t <- 1 # diffusion time
    K <- min(length(eigenvalues), ncol(eigenvectors)) # use available components
    # Precompute squared eigenvalues scaled by diffusion time
    lambda_t <- (eigenvalues[seq_len(K)])^(2 * t)
    # Initialize an empty matrix for the diffusion distances
    n <- nrow(eigenvectors)
    index_pairs <- which(upper.tri(matrix(1, n, n)), arr.ind = TRUE)
    # Calculate the pairwise diffusion distance
    distances <- apply(index_pairs, 1, .calDif,
        lambda_t = lambda_t,
        eigenvector = eigenvectors, K = K
    )
    D_diffusion <- matrix(0, n, n)
    D_diffusion[upper.tri(D_diffusion)] <- distances
    D_diffusion <- D_diffusion + t(D_diffusion)
    # `D_diffusion` is now the diffusion distance matrix
    Dif <- D_diffusion[, waypoints]
    D_ravel <- as.vector(Dif)
    sdv <- sd(D_ravel) * 1.06 * (length(D_ravel)^(-1 / 5))
    W <- exp(-0.5 * ((Dif / sdv)^2))
    W <- W / apply(W, 1, sum)
    prob <- W %*% probabilities
    if (verbose) message("Complete.")
    return(prob)
}


#' function help to calculate the diffusion distance
#'
#' @param idx integer the index of the calculated value
#' @param eigenvector numeric vector, the eigenvector from diffusion map
#' @param lambda_t eigenvalues to the power of t(diffusion time)
#' @param K The number of the eigenvectors to be used in calculation
#' @keywords internal
#' @return updated diffusion distance matrix after one iteration
.calDif <- function(idx, eigenvector, lambda_t, K) {
    i <- idx[1]
    j <- idx[2]
    diff <- eigenvector[i, seq_len(K)] - eigenvector[j, seq_len(K)]
    sqrt(sum(lambda_t * diff^2))
}
