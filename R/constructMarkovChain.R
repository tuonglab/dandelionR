#' .constructMarkovChain
#'
#' Markov chain construction
#' @param wp_data Multi scale data of the waypoints
#' @param knn. Number of nearest neighbors for graph construction
#' @param pseudotime pseudotime ordering of cells
#' @param waypoints integer vector, index of selected waypoint used to
#' @param vb whether to print messages
#' @param use_RANN parameter to make user choose
#' whether to use RANN to construct Markov chain,
#' or keep using bluster
#'
#' @keywords internal
#' @importFrom bluster makeKNNGraph
#' @importFrom igraph ends E<- E as_adjacency_matrix
#' @importFrom Matrix sparseMatrix summary
#' @importFrom purrr pmap
#' @importFrom stats dist
#' @return transition matrix of the markov chain
.constructMarkovChain <- function(wp_data, knn.,
                                  pseudotime, waypoints, vb, use_RANN) {
    if (vb) message("Markov chain construction...")
    pseudotime <- pseudotime[waypoints]
    KNNind <- if (use_RANN) {
        .RANNinx(wp_data, knn.)
    } else {
        .KNNind(wp_data, knn.)
    }
    KNN <- KNNind$KNN
    ind <- KNNind$ind
    ## select Standard deviation allowing for 'back' edges
    adaptive.k <- min(c(floor(knn. / 3) - 1, 30))
    dist_ <- lapply(ind, function(y) KNN@x[y])
    dist_sort <- lapply(dist_, sort, decreasing = TRUE)
    adaptive.std <- vapply(dist_sort, "[", adaptive.k, FUN.VALUE = double(1))
    # Directed graph construction pseudotime position of all the
    # neighbors
    traj_nbrs <- lapply(ind, function(x) pseudotime[x])
    ## Remove edges that move backwards in pseudotime except for
    ## edges that are within the computed standard deviation
    rem_edges <- pmap(list(
        .x = traj_nbrs, .y = (pseudotime - adaptive.std), .z = ind
    ), function(.x, .y, .z) {
        .z[.x < .y]
    })
    KNN <- Reduce(function(Knn, i) {
        return(.removeEdge(Knn, i, rem_edges = rem_edges))
    }, seq_len(length(waypoints)), init = KNN)
    # determine the indices and update adjacency matrix
    # cell_mapping <- seq_len(length(waypoints))
    # seems cell_mapping is not used, if there is an error them edit back
    ids <- summary(KNN)
    # anisotropic Diffusion Kernel
    aff <- exp(-(ids$x^2) / (adaptive.std[ids$i]^2) * 0.5 -
        (ids$x^2) / (adaptive.std[ids$j]^2) * 0.5)
    W <- sparseMatrix(
        i = ids$i, j = ids$j, x = aff,
        dims = dim(KNN), giveCsparse = TRUE
    )
    # Transition matrix
    D <- apply(W, 1, sum)
    ids <- summary(W)
    T_ <- sparseMatrix(
        i = ids$i, j = ids$j, x = ids$x / D[ids$i],
        dims = dim(KNN), giveCsparse = TRUE
    )
    return(T_)
}

#' Calculate the weighted adjacency matrix of knn graph and its index
#' @param wp_data Multi scale data of the waypoints
#' @param knn. Number of nearest neighbors for graph construction
#' @keywords internal
#' @importFrom bluster makeKNNGraph
#' @importFrom igraph ends E<- E as_adjacency_matrix
#' @importFrom stats dist
#' @return a list containing the weight adjacent matrix and index
.KNNind <- function(wp_data, knn.) {
    nbrs <- makeKNNGraph(wp_data, k = knn.)
    ## calculate distance of each edge
    distance_m <- as.matrix(dist(wp_data))
    edge_indices <- ends(nbrs, E(nbrs))
    weights <- vapply(seq_len(nrow(edge_indices)), function(i) {
        ### Get the nodes connected by each edge
        node1 <- as.numeric(edge_indices[i, 1])
        node2 <- as.numeric(edge_indices[i, 2])
        ### Use the distance matrix to get the distance between the
        ### two nodes
        distance_m[node1, node2]
    }, numeric(1))
    E(nbrs)$weight <- weights
    ## generate weighted adjacent matrix of this knn graph
    KNN <- as_adjacency_matrix(nbrs, attr = "weight")
    ## generate the index of each neighbor
    knn_summary <- summary(KNN)
    ind <- split(knn_summary$j, knn_summary$i)
    return(list(KNN = KNN, ind = ind))
}

#' function used in Reduce to remove KNN's backward edges except for
#'  edges that are within the computed standard deviation
#'
#' @param Knn weight KNN adjacent matrix
#' @param i the iteration number
#' @param rem_edges the edges that need to be removes
#' @keywords internal
#' @return an updated matrix after one round of iteration
.removeEdge <- function(Knn, i, rem_edges) {
    if (length(Knn[i, rem_edges[[i]]])) {
        Knn[i, rem_edges[[i]]] <- 0
    }
    return(Knn)
}

#' Calculate the weight adjacent matricks of knn graph and its index using RANN
#' @param wp_data Multi scale data of the waypoints
#' @param knn. Number of nearest neighbors for graph construction
#' @keywords internal
#' @importFrom RANN nn2
#' @importFrom Matrix sparseMatrix
#' @importFrom stats dist
#' @return a list containing the weight adjacent matrix and index
.RANNinx <- function(wp_data, knn.) {
    RAKNN <- nn2(wp_data, k = knn.)
    row_indices <- rep(seq_len(nrow(RAKNN[["nn.dists"]])),
        times = ncol(RAKNN[["nn.dists"]])
    )
    col_indices <- as.vector(RAKNN[["nn.idx"]])
    x_values <- as.vector(RAKNN[["nn.dists"]])
    KNN <- sparseMatrix(
        i = row_indices,
        j = col_indices,
        x = x_values
    )
    ## generate the index of each neighbor
    ind <- split(
        as.vector(RAKNN[["nn.idx"]]),
        rep(1:nrow(RAKNN[["nn.idx"]]),
            times = ncol(RAKNN[["nn.idx"]])
        )
    )
    return(list(KNN = KNN, ind = ind))
}
