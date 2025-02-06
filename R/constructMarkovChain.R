#' .constructMarkovChain
#'
#' Markov chain construction
#' @param wp_data Multi scale data of the waypoints
#' @param knn. Number of nearest neighbors for graph construction
#' @param pseudotime pseudotime ordering of cells
#' @param waypoints integer vector, index of selected waypoint used to
#' @param verbose whether to print messages
#' construct markov chain
#' @keywords internal
#' @importFrom bluster makeKNNGraph
#' @importFrom igraph ends E<- E as_adjacency_matrix
#' @importFrom Matrix sparseMatrix summary
#' @importFrom purrr pmap
#' @importFrom stats dist
#' @return transition matrix of the markov chain
.constructMarkovChain <- function(wp_data, knn., pseudotime, waypoints, verbose = TRUE) {
    if (verbose) {
        message("Markov chain construction...")
    }
    pseudotime <- pseudotime[waypoints]
    KNNind <- .KNNind(wp_data, knn.)
    KNN <- KNNind$KNN
    ind <- KNNind$ind
    idx_seq <- KNNind$idx_seq
    ## select Standard deviation allowing for 'back' edges
    adaptive.k <- min(c(floor(knn. / 3) - 1, 30))
    dist_ <- lapply(idx_seq, function(y) {
        KNN@x[y]
    })
    dist_sort <- lapply(dist_, sort, decreasing = TRUE)
    adaptive.std <- vapply(dist_sort, "[", adaptive.k, FUN.VALUE = double(1))
    # Directed graph construction pseudotime position of all the
    # neighbors
    traj_nbrs <- lapply(ind, function(x) {
        pseudotime[x]
    })
    ## Remove edges that move backwards in pseudotime except for
    ## edges that are within the computed standard deviation
    rem_edges <- pmap(list(
        .x = traj_nbrs, .y = (pseudotime - adaptive.std),
        .z = ind
    ), function(.x, .y, .z) {
        .z[.x < .y]
    })
    for (i in seq_len(length(waypoints))) {
        if (length(KNN[i, rem_edges[[i]]])) {
            KNN[i, rem_edges[[i]]] <- 0
        }
    }
    # determine the indice and update adjacency matrix
    # cell_mapping <- seq_len(length(waypoints))
    # seems cell_mapping is not used, if there is an error them edit back
    requireNamespace("Matrix")
    ids <- Matrix::summary(KNN)
    # anisotropic Diffusion Kernel
    aff <- exp(-(ids$x^2) / (adaptive.std[ids$i]^2)
        * 0.5 - (ids$x^2) / (adaptive.std[ids$j]^2) * 0.5)
    W <- Matrix::sparseMatrix(
        i = ids$i, j = ids$j, x = aff,
        dims = dim(KNN), giveCsparse = TRUE
    )
    # Transition matrix
    D <- apply(W, 1, sum)
    ids <- Matrix::summary(W)
    T_ <- sparseMatrix(
        i = ids$i, j = ids$j, x = ids$x / D[ids$i],
        dims = dim(KNN), giveCsparse = TRUE
    )
    return(T_)
}

#' Calculate the weight adjacent matricks of knn graph and its index
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
    edge_indices <- ends(nbrs, igraph::E(nbrs))
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
    idx <- KNN@p
    idx_seq <- mapply(seq, (idx + 1)[-length(idx)], idx[-1])
    ind <- lapply(idx_seq, function(x) {
        KNN@i[x] + 1
    })
    return(list(KNN = KNN, ind = ind, idx_seq = idx_seq))
}
