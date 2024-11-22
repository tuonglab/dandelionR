#' project_pseudotime_to_cell
#'
#' Function to project pseudotime & branch probabilities from pb.milo (pseudobulk) to adata (cell).
#' @param adata SingleCellExperiment object
#' @param pb_adata pseudobulk data
#' @param term_states vector of terminal states with branch_probabilities to be transferred
#' @param suffix suffix to be added after the added column names, default ''
#' @return subset of adata where cell that do not belong to any neighbourhood are removed and projected pseudotime information stored colData
#' @import miloR
#' @import SingleCellExperiment
#' 
#' @export
project_pseudotime_to_cell <- function(adata, pb_adata, term_states, suffix = "") {
  nhood <- nhoods(pb_adata)  # peudobulk x cells
  # leave out cells that do not blongs to any neighbourhood
  nhoodsum <- apply(nhoods(pb_adata), 2, sum)
  cdata <- adata[, nhoodsum > 0]
  message(sprintf("%d number of cells removed due to not belonging to any neighbourhood", sum(nhoodsum == 0)))
  # for each cell pesudotime_mean is the average of the pseudobulks the cell is
  # in, weighted by 1/ neighbourhood size
  nhoods_cdata <- nhood[, nhoodsum > 0]
  nhoods_cdata_norm <- nhoods_cdata/apply(nhoods_cdata, 1, sum)
  col_list <- c(paste0("pseudotime", suffix), paste0(names(term_states), suffix))
  new_col <- vapply(colData(pb_adata)[, col_list]@listData, function(x, y) {
    list(as.vector(x %*% y/apply(y, 2, sum)))
  }, y = nhoods_cdata_norm, FUN.VALUE = list(double()))
  colData(cdata) <- cbind(colData(cdata), new_col)
  cdata
}