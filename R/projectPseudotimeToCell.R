#' Project Pseudotime and Branch Probabilities to Single Cells
#'
#' This function projects pseudotime and branch probabilities from pseudobulk (`pb.milo`)
#' data to single-cell resolution (`milo`). The results are stored in the `colData` of
#' the `milo` object.
#'
#' @param milo A `SingleCellExperiment` or `Milo` object. Represents single-cell data where
#' pseudotime and branch probabilities will be projected.
#' @param pb_milo A pseudobulk `Milo` object. Contains aggregated branch probabilities and
#' pseudotime information to be transferred to single cells.
#' @param term_states A named vector of terminal states, with branch probabilities to be
#' transferred. The names should correspond to branches of interest.
#' @param pseudotime_key Character. The column name in `colData` of `pb_milo` that contains
#' the  pseudotime information which was used in the `markovProbability` function. Default
#' is `"pseudotime"`.
#' @param suffix Character. A suffix to be added to the new column names in `colData`.
#' Default is an empty string (`''`).
#' @param verbose Boolean, whether to print messages/warnings.
#'
#' @examples
#' data(sce_vdj)
#' sce_vdj <- setupVdjPseudobulk(sce_vdj,
#'     already.productive = FALSE,
#'     allowed_chain_status = c("Single pair", "Extra pair")
#' )
#' # Build Milo Object
#' set.seed(100)
#' milo_object <- miloR::Milo(sce_vdj)
#' milo_object <- miloR::buildGraph(milo_object, k = 50, d = 20, reduced.dim = "X_scvi")
#' milo_object <- miloR::makeNhoods(milo_object, reduced_dims = "X_scvi", d = 20)
#'
#' # Construct Pseudobulked VDJ Feature Space
#' pb.milo <- vdjPseudobulk(milo_object, col_to_take = "anno_lvl_2_final_clean")
#' pb.milo <- scater::runPCA(pb.milo, assay.type = "Feature_space")
#'
#' # Define root and branch tips
#' pca <- t(as.matrix(SingleCellExperiment::reducedDim(pb.milo, type = "PCA")))
#' branch.tips <- c(189, 198) # which.min(pca[, 2]) and which.max(pca[, 2])
#' names(branch.tips) <- c("CD8+T", "CD4+T")
#' root <- 177 # which.max(pca[, 1])
#'
#' # Construct Diffusion Map
#' dm <- destiny::DiffusionMap(t(pca), n_pcs = 50, n_eigs = 10)
#' dif.pse <- destiny::DPT(dm, tips = c(root, branch.tips), w_width = 0.1)
#'
#' # Markov Chain Construction
#' pb.milo <- markovProbability(
#'     milo = pb.milo,
#'     diffusionmap = dm,
#'     terminal_state = branch.tips,
#'     diffusiontime = dif.pse[[paste0("DPT", root)]],
#'     root_cell = root,
#'     pseudotime_key = "pseudotime"
#' )
#'
#' # Project Pseudobulk Data
#' projected_milo <- projectPseudotimeToCell(
#'     milo_object,
#'     pb.milo,
#'     branch.tips,
#'     pseudotime_key = "pseudotime"
#' )
#'
#' @return subset of milo or SingleCellExperiment object where cell that do not belong to any neighbourhood are removed and projected pseudotime information stored colData
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData<-
#' @importFrom S4Vectors metadata
#' @export
projectPseudotimeToCell <- function(milo, pb_milo, term_states = NULL, pseudotime_key = "pseudotime", suffix = "", verbose = TRUE) {
    if (is.null(term_states)) {
        if (is.null(metadata(pb_milo)$branch.tips)) # nocov start
            {
                abort("Parameter Error: Please provide term_state, which should align with parameter terminal_state in function markovProbability")
            } # nocov end
        else {
            term_states <- metadata(pb_milo)$branch.tips
        }
    }
    nhood <- miloR::nhoods(pb_milo) # peudobulk x cells
    # leave out cells that do not blongs to any neighbourhood
    nhoodsum <- apply(miloR::nhoods(pb_milo), 2, sum)
    cdata <- milo[, nhoodsum > 0]
    if (verbose) {
        message(sprintf(
            "%d number of cells removed due to not belonging to any neighbourhood",
            sum(nhoodsum == 0)
        ))
    }
    # for each cell pesudotime_mean is the average of the pseudobulks the cell
    # is in, weighted by 1/ neighbourhood size
    nhoods_cdata <- nhood[, nhoodsum > 0]
    nhoods_cdata_norm <- nhoods_cdata / apply(nhoods_cdata, 1, sum)
    col_list <- c(paste0(pseudotime_key, suffix), paste0(names(term_states), suffix))
    new_col <- vapply(colData(pb_milo)[, col_list]@listData, function(x, y) {
        list(as.vector(x %*% y / apply(y, 2, sum)))
    }, y = nhoods_cdata_norm, FUN.VALUE = list(double()))
    colData(cdata) <- cbind(colData(cdata), new_col)
    return(cdata)
}
