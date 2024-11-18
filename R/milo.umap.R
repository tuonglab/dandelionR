#' milo.umap
#' 
#' conduct the umap on the knn graph in umap
#' @param milo the milo object need to conduct umap on
#' @param slot_name the slot name in reducedim where the result store
#' @param n.neighbors the number of neighboring points used
#' @param metric the choice of metric used to measure distance 
#' @return milo object with umap reduction
#' @export
milo.umap <- function(milo,slot_name = "UMAP_graph", n.neighbors = 50L, metric = "euclidean")
{
  # get the graph's adjacency matrix
  graphm <- as_adjacency_matrix(milo@graph[[1]])
  # inherit the names of each row
  rownames(graphm) <- rownames(colData(milo))
  colnames(graphm)<- rownames(colData(milo))
  
  # conduct umap
  umap <- RunUMAP(graphm, n.neighbors = n.neighbors, metric = metric)
  reducedDim(milo,"UMAP_graph") <- umap@cell.embeddings
  milo
}