#' milo_umap
#'
#' use function RunUMAP from Seurat to conduct the umap on the adjacency matrix of knn graph in milo object
#' @param milo the milo object with knn graph that needed to conduct umap on.
#' @param slot_name character, with default 'UMAP_knngraph'.
#'  - The slot name in reduceDim where the result store
#' @param n.neighbors integer, with default 50L.
#'  - the number of neighboring points used
#'  - parameter of RunUMAP, checking its document for further detail
#' @param metric character, with default 'euclidean'
#'  - the choice of metric used to measure distance
#'  - parameter of RunUMAP, checking its document for further detail
#' @return milo object with umap reduction
#' @import SingleCellExperiment
#' @examples
#' 
#' # load denpendency 
#' library(miloR)
#' 
#' # load example data
#' data(exmilo)
#' 
#' # Construct UMAP
#' exmilo <- milo_umap(exmilo, n.neighbors = 10L, metric = "euclidean")
#' 
#' # visualize the result
#' scater::plotUMAP(exmlio,dimred = "UMAP_knngraph")
#' 
#' @export
milo_umap <- function(milo, slot_name = "UMAP_knngraph", n.neighbors = 50L, metric = "euclidean") {
    requireNamespace("miloR")
    requireNamespace("igraph")
    # get the graph's adjacency matrix
    graphm <- igraph::as_adjacency_matrix(miloR::graph(milo))
    # inherit the names of each row
    rownames(graphm) <- rownames(colData(milo))
    colnames(graphm) <- rownames(colData(milo))
    # conduct umap
    requireNamespace("Seurat")
    umap <- Seurat::RunUMAP(graphm, n.neighbors = n.neighbors, metric = metric)
    reducedDim(milo, slot_name) <- umap@cell.embeddings
    milo
}
