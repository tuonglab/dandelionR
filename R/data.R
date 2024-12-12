#' Example Dataset for V(D)J Analysis
#'
#' The `sce_vdj` object is a down-sampled demo dataset derived from Suo et al., 2024, *Nature Biotechnology*.
#' This dataset is used in vignettes to demonstrate workflows for V(D)J analysis. For details, see the original
#' publication at \url{https://www.nature.com/articles/s41587-023-01734-7}.
#'
#' @docType data
#' @usage data(sce_vdj)
#' @format A `SingleCellExperiment` object with the following slots:
#' \describe{
#'   \item{\code{colData}}{
#'     A `DataFrame` containing metadata about each sample, corresponding to `obs` in AnnData (Python).
#'     The following columns are relevant for vignette usage:
#'     \describe{
#'       \item{\code{productive_(mode)_VDJ}, \code{productive_(mode)_VJ}}{
#'         Factors indicating whether the heavy or light chain is productive. \code{mode} refers to the
#'         extraction mode for V(D)J genes and can be one of:
#'         \itemize{
#'           \item{\code{'abT'}}: TCRαβ
#'           \item{\code{'gdT'}}: TCRγδ
#'           \item{\code{'B'}}: BCR
#'         }
#'       }
#'       \item{\code{Gene segment fields}}{
#'         Gene segment annotations with column names in the format \code{(v/d/j)_call_(mode)_(VDJ/VJ)}.
#'         Examples include:
#'         \itemize{
#'           \item{\code{v_call_abT_VDJ}}: V gene for TCRαβ VDJ recombination
#'           \item{\code{d_call_abT_VJ}}: D gene for TCRαβ VJ recombination
#'         }
#'       }
#'       \item{\code{chain_status}}{
#'         A factor describing the receptor chain's status.
#'       }
#'       \item{\code{anno_lvl_2_final_clean}}{
#'         Cell type annotations.
#'       }
#'     }
#'   }
#'   \item{\code{int_colData}}{
#'     A `DataFrame` containing additional assay metadata important for further analysis. Includes:
#'     \itemize{
#'       \item{\code{X_scvi}}: A dimensionality reduction matrix from the scVI model.
#'       \item{\code{UMAP}}: A UMAP reduction matrix.
#'     }
#'   }
#' }
#'
#' @source Suo et al., 2024, *Nature Biotechnology*. \url{https://www.nature.com/articles/s41587-023-01734-7}.
#' @examples
#' data(sce_vdj)
"sce_vdj"

#' Example SCE Dataset that does not contain V(D)J information
#'
#' The `demo_sce` object is a down-sampled demo dataset derived from Suo et al., 2024, *Nature Biotechnology*.
#' This dataset is used in vignettes to demonstrate workflows for V(D)J analysis. For details, see the original
#' publication at \url{https://www.nature.com/articles/s41587-023-01734-7}. The original Lymphoid cells data in
#' h5ad format is available at \url{https://developmental.cellatlas.io/fetal-immune}.
#'
#' @docType data
#' @usage data(demo_sce)
#' @format A `SingleCellExperiment` object with the following slots:
#' \describe{
#'   \item{\code{colData}}{
#'     A minimall `DataFrame` containing metadata about each sample, corresponding to `obs` in AnnData (Python).
#'     The following columns are relevant for vignette usage:
#'     \describe{
#'       \item{\code{anno_lvl_2_final_clean}}{
#'         Cell type annotations.
#'       }
#'     }
#'   }
#'   \item{\code{int_colData}}{
#'     A `DataFrame` containing additional assay metadata important for further analysis. Includes:
#'     \itemize{
#'       \item{\code{X_scvi}}: A dimensionality reduction matrix from the scVI model.
#'       \item{\code{UMAP}}: A UMAP reduction matrix.
#'     }
#'   }
#' }
#'
#' @source Suo et al., 2024, *Nature Biotechnology*. \url{https://www.nature.com/articles/s41587-023-01734-7}.
#' @examples
#' data(demo_sce)
"demo_sce"

#' Example AIRR Dataset for V(D)J Analysis
#'
#' The `demo_airr` object is a list of AIRR data frames from a down-sampled demo dataset derived from Suo et al.,
#' 2024, *Nature Biotechnology*. This dataset is used in vignettes to demonstrate workflows for V(D)J analysis.
#' For details, see the original publication at \url{https://www.nature.com/articles/s41587-023-01734-7}. The
#' original files are available at \url{https://github.com/zktuong/dandelion-demo-files}.
#'
#' @docType data
#' @usage data(demo_airr)
#' @format A `SingleCellExperiment` object with the following slots:
#' \describe{
#'   \item{\code{list}}{
#'     List of `DataFrames` containing the standardised AIRR data for each sample.
#'     For information of AIRR rearrangements, see the AIRR Community standards at
#'     \url{https://docs.airr-community.org/en/stable/datarep/rearrangements.html}.
#'   }
#' }
#'
#' @source Suo et al., 2024, *Nature Biotechnology*. \url{https://www.nature.com/articles/s41587-023-01734-7}.
#' @examples
#' data(demo_airr)
"demo_airr"
