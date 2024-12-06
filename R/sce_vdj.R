#' demo data used in vignettes
#'
#' sce_vdj - A down-sampled demo data from Suo et al 2024 Nature Biotechnology. See \url{https://www.nature.com/articles/s41587-023-01734-7}
#' @docType data
#' @usage data(sce_vdj)
#' @format A SingleCellExperiment object with the following slots filled
#' \describe{
#'  \item{\code{colData}}{
#'    A \code{DataFrame} containing metadata about each sample, corresponding to obs in AnnData for python.\cr
#'    Below only describe the column used in the vignette:
#'    \describe{
#'      \item{\code{productive_(mode)_VDJ}, \code{productive_(mode)_VJ}}{
#'        factor containing the information about whether the heavy/light chain is productive
#'        mode: mode for extraction of the V(D)J genes, could be 'abT'(TCRαβ), 'gdT'(TCRγδ) or 'B'(BCR)
#'      }
#'      \item{Gene segment fields}{
#'        with a pattern like \code{(v/d/j)_call_(mode)_(VDJ/VJ)}\cr For example:
#'        \itemize{
#'           \item{\code{v_call_abT_VDJ}}: V gene for TCRαβ VDJ recombination
#'           \item{\code{d_call_abT_VJ}}: D gene for TCRαβ VJ recombination
#'        }
#'      }
#'      \item{\code{chain_status}}{the status of the receptor chain}
#'      \item{\code{anno_lvl_2_final_clean}}{cell type annotation}
#'    }
#'  }
#'  \item{\code{int_colData}}{
#'    A \code{DataFrame} contains metadata of sample, mostly additional assay that is important for further analysis.\cr
#'    This int_colData only has the reduction matrix
#'    \itemize{
#'      \item{\code{X_scvi}}: result of scVI model reduction
#'      \item{\code{UMAP}}: result of umap reduction
#'    }
#'  }
#' }
#' @source \url{https://www.nature.com/articles/s41587-023-01734-7}
#' @examples
#' data(sce_vdj)
"sce_vdj"
