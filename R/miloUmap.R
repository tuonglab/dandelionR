#' Perform UMAP on the Adjacency Matrix of a Milo Object
#'
#' This function uses `uwot::umap` to perform UMAP dimensionality
#' reduction on the adjacency matrix of the KNN graph in a Milo object.
#'
#' @param milo the milo object with knn graph that needed to conduct umap on.
#' @param slot_name character, with default 'UMAP_knngraph'.
#'  - The slot name in reduceDim where the result store
#' @param n_neighbors integer, with default 50L.
#'  - the size of local neighborhood (in terms of number of
#'   neighboring sample points) used for manifold approximation.
#'  - Here, the goal is to create large enough neighborhoods to capture
#'   the local manifold structure to allow for hypersampling.
#' @param metric character, with default 'euclidean'
#'  - the choice of metric used to measure distance to find nearest neighbors.
#'   Default is 'euclidean'.
#' @param min_dist numeric, with default 0.3
#'  - the minimum distance between points in the low dimensional space
#' @param use_graph Logical, default TRUE.
#'   - Whether to run UMAP on the graph adjacency matrix (TRUE) as in Dandelion,
#'     or directly on the latent space (FALSE) for faster performance.
#' @param ... other parameters passed to uwot::umap
#' @examples
#' data(sce_vdj)
#' # downsample to just 1000 cells
#' sce_vdj <- sce_vdj[, 1:1000]
#' sce_vdj <- setupVdjPseudobulk(sce_vdj,
#'     already.productive = FALSE,
#'     allowed_chain_status = c("Single pair", "Extra pair")
#' )
#' # Build Milo Object
#' milo_object <- miloR::Milo(sce_vdj)
#' milo_object <- miloR::buildGraph(milo_object,
#'     k = 50, d = 20,
#'     reduced.dim = "X_scvi"
#' )
#' milo_object <- miloR::makeNhoods(milo_object,
#'     reduced_dims = "X_scvi", d = 20
#' )
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
miloUmap <- function(
    milo, slot_name = "UMAP_knngraph", n_neighbors = 50L, metric = "euclidean",
    min_dist = 0.3, use_graph = TRUE,...) {
    if(use_graph)
    {
      # get the graph's adjacency matrix
      graphm <- as_adjacency_matrix(graph(milo), sparse = TRUE)
      # inherit the names of each row
      rownames(graphm) <- rownames(colData(milo))
      colnames(graphm) <- rownames(colData(milo))
    }
    else
    {
      graph <- reducedDim(milo_object, "X_scvi")
    }
    # conduct umap
    pos <- umap(graphm,
        n_neighbors = n_neighbors, metric = metric, min_dist = min_dist,
        ...
    )
    reducedDim(milo, slot_name) <- pos
    return(milo)
}
