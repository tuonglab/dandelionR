#' demo data used in vignettes 
#'
#' sce_vdj - A down-sampled demo data from Suo et al 2024 Nature Biotechnology. See \url{https://www.nature.com/articles/s41587-023-01734-7}
#' @docType data
#' @usage data(sce_vdj)
#' @format A SingleCellExperiment object with the following slots filled
#' \describe{
#'   \item{assays}{
#'   \itemize{Currently only contains 'counts' and 'logcounts'
#'   \item{X - Normalized and log1p transformed expression data}
#' }
#'   \item{colData}{Cell level metadata}
#' }
#' @source \url{https://www.nature.com/articles/s41587-023-01734-7}
#' @examples
#' data(sce_vdj)
"sce_vdj"