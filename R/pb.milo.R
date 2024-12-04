#' demo data used in vignettes
#'
#' pb.milo - A milo object used to storepseudobulk vdj feature space(a matrix with a shape vdj x pseudobulk)
#' @docType data
#' @usage data(pb.milo)
#' @format A Milo object with the following slots filled
#' \describe{
#'  \item{\code{Feature_space}}{VDJ feature space, a vdj x pseudobulk matrix store each vdj' usage in each pseudobulk}
#'  \item{\code{nhoods}}{a pseudobulk x cell sparse matrix store the information of whether the cell is contained in this pseudobulk, 1 ture, 0 false.}
#'  \item{\code{colData}}{
#'    A \code{DataFrame} containing metadata about each sample, corresponding to obs in AnnData for python.
#'    \describe{
#'      \item{\code{anno_lvl_2_final_clean}}{most common cell type in this pseudobulk}
#'      \item{\code{ano_lvl_2_final_clean}}{the exact fraction of the most common cell type}
#'      \item{\code{cell_count}}{the number of cell in each pseudobulk}
#'    }
#'  }
#' }
#' @examples
#' data(pb.milo)
"pb.milo"