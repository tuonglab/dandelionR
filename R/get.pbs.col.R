#' .get.pbs.col
#'
#' Helper function to create the new pseudobulk object's coldata.
#' @param pbs dgeMatrix, cell x pseudobulk binary matrix
#' @param col_to_take character vector, names of colData of milo that need to be processed  
#' @param milo Milo or SingleCellExperiment object
#' @return pbs_col, a DataFrame which will be passed to the new SingleCellExperiment object as colData of vdj x pseudobulk assays
#' @import Matrix
#' @import S4Vectors
.get.pbs.col <- function(pbs, col_to_take, milo) {
  # prepare per-pseudobulk calls of specified metadata columns
  if (!is.null(col_to_take)) {
    pbs.col <- DataFrame()
    for (anno_col in col_to_take) {
      fa <- as.factor(colData(milo)[[anno_col]])
      fa <- data.frame(fa)
      anno.dummies <- model.matrix(~. - 1, data = fa, contrasts.arg = lapply(fa,
        contrasts, contrasts = FALSE))
      anno.count <- Matrix::t(pbs) %*% anno.dummies  # t(cell x  pseudo)   %*% (cell x element) = pseudo x element
      anno.sum <- apply(anno.count, 1, sum)
      anno.frac <- anno.count/anno.sum
      anno.frac <- DataFrame(anno.frac)
      colnames(anno.frac) <- colnames(anno.dummies)
      pbs.col[anno_col] <- colnames(anno.frac)[apply(anno.frac, 1, which.max)]
      pbs.col[paste0(anno_col, "_fraction")] <- apply(anno.frac, 1, max)
    }
  }
  # report the number of cells for each pseudobulk ensure pbs is an array so
  # that it sums into a vector that can go in easily
  pbs.col["cell_count"] <- apply(pbs, 2, sum)
  return(pbs.col)
}