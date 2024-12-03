#' differentiation_entropy
#'
#' function to compute branch probabilities using markov chain
#' @param wp_data Multi scale data of the waypoints
#' @param terminal_states integer, NULL by default. The index of the terminal state.
#' @param knn. integer, 30L by default. Number of nearest neighbors for graph construction.
#' @param pseudotime pseudotime ordering of cells
#' @param waypoints integer vector, index of selected waypoint used to construct markov chain
#' @return probabilities
#' @include construct.markov.chain.R
#' @include terminal.state.from.markov.chain.R
differentiation_probabilities <- function(wp_data, terminal_states = NULL, knn. = 30L,
    pseudotime, waypoints) {
    T_ <- .construct.markov.chain(wp_data, 30, pseudotime, waypoints)
    # identify terminal states if not specified
    if (is.null(terminal_states)) {
        terminal_states <- .terminal.state.from.markov.chain(T_, wp_data, pseudotime,
            waypoints)
    }
    abs_states_idx <- which(waypoints %in% terminal_states)
    T_[abs_states_idx, ] <- 0
    T_ <- Reduce(function(matri, x) {
        matri[x, x] <- 1
        matri
    }, x = abs_states_idx, init = T_)
    message("Computing fundamental matrix and absorption probabilities...")
    # Transition states
    trans_states_idx <- which(!(waypoints %in% terminal_states))
    # Q matrix
    Q <- T_[-abs_states_idx, -abs_states_idx]
    # Fundamental matrix
    mat <- diag(dim(Q)[[1]]) - Q
    requireNamespace("Matrix")
    requireNamespace("MASS")
    N <- tryCatch({
        Matrix::solve(mat)
    }, error = function(cnd) {
        warning("Matrix generated is singular or nearly singular; using pseudo-inverse to construct fundamental matrix")
        warning("Or you can re-run this function to reconstruct the markov chain")
        MASS::ginv(as.matrix(mat))
    })
    R <- T_[trans_states_idx, abs_states_idx]
    # absorbing probabilities:
    probabilities <- N %*% R
    probabilities@x[probabilities@x < 0] <- 0
    # add terminal states
    probabilities <- rbind(probabilities, T_[abs_states_idx, abs_states_idx])
    probabilities <- probabilities[order(c(trans_states_idx, abs_states_idx)), ]
    return(probabilities)
}

