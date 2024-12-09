#' Perform UMAP on the Adjacency Matrix of a Milo Object
#'
#' This function uses `umap` from the `uwot` package to perform UMAP dimensionality
#' reduction on the adjacency matrix of the KNN graph in a Milo object.
#'
#' @param milo A `Milo` object containing a KNN graph. UMAP will be conducted on
#' the adjacency matrix of this graph.
#' @param slot_name Character. The name of the slot in `reducedDims` where the UMAP
#' results will be stored. Default is `'UMAP_knngraph'`.
#' @param n_neighbors Integer. The size of local neighborhood (in terms of number of
#' neighboring sample points) used for manifold approximation. Default is `50L`.
#' For further details, refer to the `uwot::umap` documentation.
#' @param metric Character. The distance metric to be used in the UMAP algorithm.
#' Default is `'euclidean'`. For further details, refer to the `uwot::umap` documentation.
#' @param min_dist Float. The minimum distance between points in the low dimensional space.
#' Default is `0.3`.
#' @param ... other parameters passed to `uwot::umap`.
#' @examples
#' data(sce_vdj)
#' sce_vdj <- setupVdjPseudobulk(sce_vdj,
#'     already.productive = FALSE
#' )
#' # Build Milo Object
#' traj_milo <- miloR::Milo(sce_vdj)
#' milo_object <- miloR::buildGraph(traj_milo, k = 50, d = 20, reduced.dim = "X_scvi")
#' milo_object <- miloR::makeNhoods(milo_object, reduced_dims = "X_scvi", d = 20)
#'
#' # Construct UMAP on Milo Neighbor Graph
#' milo_object <- miloUmap(milo_object)
#'
#' @return milo object with umap reduction
#' @import SingleCellExperiment
#' @importFrom igraph as_adjacency_matrix
#' @importFrom miloR graph
#' @importFrom uwot umap
#' @export
miloUmap <- function(milo, slot_name = "UMAP_knngraph", n_neighbors = 50L, metric = "euclidean", min_dist = 0.3, ...) {
    # get the graph's adjacency matrix
    graphm <- as_adjacency_matrix(graph(milo))
    # inherit the names of each row
    rownames(graphm) <- rownames(colData(milo))
    colnames(graphm) <- rownames(colData(milo))
    # conduct umap
    pos <- umap(graphm, n_neighbors = n_neighbors, metric = metric, min_dist = min_dist, ...)
    reducedDim(milo, slot_name) <- pos
    milo
}
