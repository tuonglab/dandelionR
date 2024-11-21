#' .construct.markov.chain
#'
#' markov chain construction
#' @param wp_data Multi scale data of the waypoints
#' @param knn. Number of nearest neighbors for graph construction
#' @param pseudotime pseudotime ordering of cells
#' @param waypoints integer vector, index of selected waypoint used to construct markov chain
#' @return transition matrix of the markov chain
#' @import bluster
#' @import igraph
#' @import Matrix
#' @import purrr
.construct.markov.chain <- function(wp_data, knn., pseudotime, waypoints) {
  message("Markov chain construction...")
  pseudotime <- pseudotime[waypoints]
  # construct kNN graph
  nbrs <- bluster::makeKNNGraph(wp_data, k = knn.)
  ## calculate distance of each edge
  distance_m <- as.matrix(stats::dist(wp_data))
  E(nbrs)$weight <- sapply(E(nbrs), function(e) {
    ### Get the nodes connected by each edge
    nodes <- ends(nbrs, e)
    node1 <- as.numeric(nodes[1])
    node2 <- as.numeric(nodes[2])
    ### Use the distance matrix to get the distance between the two nodes
    distance_m[node1, node2]
  })
  ## generate weighted adjacent matrix of this knn graph
  KNN <- as_adjacency_matrix(nbrs, attr = "weight")
  ## generate the index of each neighbor
  idx <- KNN@p
  idx_seq <- mapply(seq, (idx + 1)[-length(idx)], idx[-1])
  ind <- lapply(idx_seq, function(x) {
    KNN@i[x] + 1
  })
  ## select Standard deviation allowing for 'back' edges
  adaptive.k <- min(c(floor(knn./3) - 1, 30))
  dist_ <- lapply(idx_seq, function(y) {
    KNN@x[y]
  })
  dist_sort <- lapply(dist_, sort, decreasing = TRUE)
  adaptive.std <- vapply(dist_sort, "[", adaptive.k, FUN.VALUE = double(1))
  # Directed graph construction pseudotime position of all the neighbors
  traj_nbrs <- lapply(ind, function(x) {
    pseudotime[x]
  })
  ## Remove edges that move backwards in pseudotime except for edges that are
  ## within the computed standard deviation
  rem_edges <- pmap(list(.x = traj_nbrs, .y = (pseudotime - adaptive.std), .z = ind),
    function(.x, .y, .z) {
      .z[.x < .y]
    })
  for (i in 1:length(waypoints)) {
    if (length(KNN[i, rem_edges[[i]]]))
      KNN[i, rem_edges[[i]]] <- 0
  }
  # determine the indice and update adjacency matrix
  cell_mapping <- 1:length(waypoints)
  ids <- summary(KNN)
  # anisotropic Diffusion Kernel
  aff <- exp(-(ids$x^2)/(adaptive.std[ids$i]^2) * 0.5 - (ids$x^2)/(adaptive.std[ids$j]^2) *
    0.5)
  W = sparseMatrix(i = ids$i, j = ids$j, x = aff, dims = dim(KNN), giveCsparse = TRUE)
  # Transition matrix
  D = apply(W, 1, sum)
  ids = summary(W)
  T_ = sparseMatrix(i = ids$i, j = ids$j, x = ids$x/D[ids$i], dims = dim(KNN),
    giveCsparse = TRUE)
  return(T_)
}