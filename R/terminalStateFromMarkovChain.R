#' .terminalStateFromMarkovChain
#'
#' @param Transmat Transition matrix
#' @param wp_data Multi scale data of the waypoints
#' @param pseudotime numeric vector, pseudotime of each pseudobulk
#' @param waypoints integer vector, waypoint selected to construct markov chain.
#' @keywords internal
#' @return terminal_state
.terminalStateFromMarkovChain <- function(Transmat, wp_data, pseudotime, waypoints) {
    message("No terminal state provided, indentification of terminal states....")
    # Identify terminal states dm_boudaries
    n <- min(dim(Transmat))
    ei <- eigen(t(Transmat))
    idx <- order(Re(ei$values), decreasing = TRUE)[seq_len(10)]
    vals <- Re(ei$values[idx])
    vecs <- ei$vectors[, idx]
    dm_boudaries <- unique(apply(Transmat, 2, which.max), apply(Transmat, 2, which.min))
    ranks <- abs(Re(vecs[, which.max(Re(vals)), drop = FALSE]))
    # cutoff and intersection with the boundary cells
    requireNamespace("stats")
    cutoff <- stats::qnorm(0.9999, mean = stats::median(ranks), sd = stats::median(abs(ranks -
        stats::median(ranks))))
    # connect components of cells beyond cutoff
    cells <- which(ranks > cutoff)
    # Find connected componets
    requireNamespace("igraph")
    graph <- igraph::graph_from_adjacency_matrix(Transmat[cells, cells],
        weighted = TRUE,
        mode = "directed"
    )
    components_g <- igraph::components(graph)
    requireNamespace("purrr")
    cells <- unlist(purrr::map(.x = seq_len(components_g$no), .f = ~ {
        nodes_in_component <- which(components_g$membership == .x)
        id <- which.max(pseudotime[nodes_in_component])
        nodes_in_component[id]
    }))
    # Nearest diffusion map boundaries
    terminal_states <- c()
    requireNamespace("spam")
    for (i in cells) {
        dists <- spam::nearest.dist(wp_data[dm_boudaries, ], wp_data[i, , drop = FALSE])
        terminal_states <- c(terminal_states, dm_boudaries[which.max(dists@entries)])
    }
    unique(terminal_states)
}
