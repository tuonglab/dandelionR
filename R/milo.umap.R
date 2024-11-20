#' milo.umap
#'
#' use function RunUMAP from Seurat to conduct the umap on the knn graph in milo object
#' @param milo the milo object with knn graph that needed to conduct umap on.
#' @param slot_name character, with default "UMAP_knngraph". 
#'    The slot name in reduceDim where the result store
#' @param n.neighbors integer, with default 50L.
#'    the number of neighboring points used
#' @param metric character, the choice of metric used to measure distance
#' @return milo object with umap reduction
#' @import igraph
#' @import SingleCellExperiment
#' @import Seurat
#' @export
milo.umap <- function(milo, slot_name = "UMAP_knngraph", n.neighbors = 50L, metric = "euclidean") {
  requireNamespace("miloR")
  # get the graph's adjacency matrix
  graphm <- as_adjacency_matrix(miloR::graph(milo))
  # inherit the names of each row
  rownames(graphm) <- rownames(colData(milo))
  colnames(graphm) <- rownames(colData(milo))
  # conduct umap
  umap <- RunUMAP(graphm, n.neighbors = n.neighbors, metric = metric)
  reducedDim(milo, slot_name) <- umap@cell.embeddings
  milo
}