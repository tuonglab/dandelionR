#' Determine terminal states using Markov chain if terminal states are not provided.
#'
#' @param Transmat Transition matrix
#' @param wp_data Multi scale data of the waypoints
#' @param pseudotime numeric vector, pseudotime of each pseudobulk
#' @param waypoints integer vector, waypoint selected to construct markov chain.
#' @param verbose Boolean, whether to print messages/warnings.
#' @keywords internal
#' @importFrom igraph graph_from_adjacency_matrix components
#' @importFrom purrr map
#' @importFrom stats median qnorm
#' @importFrom spam nearest.dist
#' @return terminal_state
.terminalStateFromMarkovChain <- function(Transmat, wp_data, pseudotime, waypoints, verbose = TRUE) {
    if (verbose) message("No terminal state provided, identification of terminal states....")
    # Identify terminal states dm_boundaries
    n <- min(dim(Transmat))
    ei <- eigen(t(Transmat))
    idx <- order(Re(ei$values), decreasing = TRUE)[seq_len(10)]
    vals <- Re(ei$values[idx])
    vecs <- ei$vectors[, idx]
    dm_boudaries <- unique(apply(Transmat, 2, which.max), apply(Transmat, 2, which.min))
    ranks <- abs(Re(vecs[, which.max(Re(vals)), drop = FALSE]))
    # cutoff and intersection with the boundary cells
    cutoff <- stats::qnorm(0.9999, mean = stats::median(ranks), sd = stats::median(abs(ranks -
        stats::median(ranks))))
    # connect components of cells beyond cutoff
    cells <- which(ranks > cutoff)
    # Find connected components
    graph <- graph_from_adjacency_matrix(Transmat[cells, cells],
        weighted = TRUE,
        mode = "directed"
    )
    components_g <- components(graph)
    cells <- unlist(map(.x = seq_len(components_g$no), .f = ~ {
        nodes_in_component <- which(components_g$membership == .x)
        id <- which.max(pseudotime[nodes_in_component])
        nodes_in_component[id]
    }))
    # Nearest diffusion map boundaries
    terminal_states <- c()
    terminal_states <- Reduce(function(term, i) {
        terminal_states <- .determTerminal(
            terminal_states = term, i = i,
            dm_boudaries = dm_boudaries, wp_data = wp_data
        )
        return(terminal_states)
    }, cells, init = terminal_states)
    unique(waypoints[terminal_states])
}

#' .determTerminal
#'
#' function in Reduce to provide waypoints
#' @param terminal_states integer vector to store the generated waypoint index
#' @param i iteration index
#' @param dm_boudaries index of the maxium or minium value of
#'  transition matrix per row
#' @param wp_data Multi scale data of the waypoints
#' @keywords internal
#' @return integer vector store the index of waypoints serve as terminal state
.determTerminal <- function(terminal_states, i, dm_boudaries, wp_data) {
    dists <- nearest.dist(wp_data[dm_boudaries, ], wp_data[i, , drop = FALSE])
    terminal_states <- c(terminal_states, dm_boudaries[which.max(dists@entries)])
    return(terminal_states)
}
