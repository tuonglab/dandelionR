#' Project Pseudotime and Branch Probabilities to Single Cells
#'
#' This function projects pseudotime and branch probabilities from pseudobulk
#' data to single-cell resolution (`milo`). The results are stored in the
#' `colData` of the `milo` object.
#'
#' @param milo A `SingleCellExperiment` or `Milo` object. Represents single-cell
#' data where pseudotime and branch probabilities will be projected.
#' @param pb_milo A pseudobulk `Milo` object. Contains aggregated branch
#' probabilities and pseudotime information to be transferred to single cells.
#' @param value_key Character. The column name in `colData` of `pb_milo`
#' that contains the value that is needed to be projected back. Default is `NULL`.
#' @param suffix Character. A suffix to be added to the new column names in
#' `colData`. Default is an empty string (`''`).
#' @param verbose Boolean, whether to print messages/warnings.
#'
#' @examples
#' data(sce_vdj)
#' # downsample to first 2000 cells
#' sce_vdj <- sce_vdj[, 1:2000]
#' sce_vdj <- setupVdjPseudobulk(sce_vdj,
#'     already.productive = FALSE,
#'     allowed_chain_status = c("Single pair", "Extra pair")
#' )
#' # Build Milo Object
#' set.seed(100)
#' milo_object <- miloR::Milo(sce_vdj)
#' milo_object <- miloR::buildGraph(milo_object,
#'     k = 50, d = 20,
#'     reduced.dim = "X_scvi"
#' )
#' milo_object <- miloR::makeNhoods(milo_object,
#'     reduced_dims = "X_scvi",
#'     d = 20
#' )
#'
#' # Construct Pseudobulked VDJ Feature Space
#' pb.milo <- vdjPseudobulk(milo_object, col_to_take = "anno_lvl_2_final_clean")
#' pb.milo <- scater::runPCA(pb.milo, assay.type = "Feature_space")
#'
#' # Define root and branch tips
#' pca <- t(as.matrix(SingleCellExperiment::reducedDim(pb.milo, type = "PCA")))
#' branch.tips <- c(which.min(pca[, 2]), which.max(pca[, 2]))
#' names(branch.tips) <- c("CD8+T", "CD4+T")
#' root <- which.min(pca[, 1])
#'
#' # Construct Diffusion Map
#' dm <- destiny::DiffusionMap(t(pca), n_pcs = 10, n_eigs = 5)
#' dif.pse <- destiny::DPT(dm, tips = c(root, branch.tips), w_width = 0.1)
#'
#' # Markov Chain Construction
#' pb.milo <- markovProbability(
#'     milo = pb.milo,
#'     diffusionmap = dm,
#'     diffusiontime = dif.pse[[paste0("DPT", root)]],
#'     terminal_state = branch.tips,
#'     root_cell = root,
#'     pseudotime_key = "pseudotime"
#' )
#' # Project Pseudobulk Data
#' projected_milo <- projectPseudotimeToCell(
#'     milo_object,
#'     pb.milo,
#'     value_key = c("pseudotime", "CD8+T", "CD4+T")
#' )
#'
#' @return subset of milo or SingleCellExperiment object where cell that do not
#' belong to any neighbourhood are removed and projected pseudotime information
#' stored colData
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData<-
#' @importFrom S4Vectors metadata
#' @export
projectPseudotimeToCell <- function(
    milo, pb_milo, value_key = NULL, suffix = "", verbose = TRUE) {
    if (is.null(value_key)) {
        abort("Please specify the column name(s) of the value(s) to be projected back.")
    } else if (!all(value_key %in% colnames(colData(pb_milo)))) {
        missing_value <- value_key[!(value_key %in% colnames(colData(pb_milo)))]
        abort(sprintf(
            "value %s do(es) not exist in pseudobulk `Milo` object.",
            paste(missing_value, collapse = ", ")
        ))
    }
    nhood <- miloR::nhoods(pb_milo) # peudobulk x cells
    # leave out cells that do not blongs to any neighbourhood
    nhoodsum <- apply(miloR::nhoods(pb_milo), 2, sum)
    cdata <- milo[, nhoodsum > 0]
    if (verbose) {
        message(sprintf(
            c(
                "%d number of cells removed due to not belonging to any",
                " neighbourhood"
            ),
            sum(nhoodsum == 0)
        ))
    }
    # for each cell pesudotime_mean is the average of the pseudobulks the cell
    # is in, weighted by 1/ neighbourhood size
    nhoods_cdata <- nhood[, nhoodsum > 0]
    nhoods_cdata_norm <- nhoods_cdata / apply(nhoods_cdata, 1, sum)
    col_list <- value_key
    new_col <- mapply(
        function(x, value_name) {
            project_single_value(x, nhoods_cdata_norm, value_name, verbose)
        },
        colData(pb_milo)[, col_list]@listData,
        col_list,
        SIMPLIFY = FALSE
    )
    colData(cdata) <- cbind(colData(cdata), new_col)
    return(cdata)
}

#' Function to project pseudobulk-level values to single-cell level
#'
#' @param x Numeric vector, pseudobulk-level value to be projected.
#' @param y Matrix (pseudobulk x cell), used to project x back to cell level.
#' @param value_name Character, name of the value being projected.
#' @param verbose Boolean, whether to print messages/warnings.
#' @return Numeric vector of projected values at cell level.
project_single_value <- function(x, y, value_name, verbose = TRUE) {
    na_idx <- is.na(x)
    if (any(na_idx)) { # nocov
        if (all(na_idx)) {
            abort(paste("All values for", value_name, "are NA."))
        }
        # exclude pseudobulk with NA value
        x <- x[!na_idx]
        modi_y <- y[!na_idx, , drop = FALSE]
        col_sums <- Matrix::colSums(modi_y)
        zero_cols <- sum(col_sums == 0)
        if (verbose) {
            message(sprintf(
                "%d cells do not have projected values for %s because their neighbour(s) value is NA",
                zero_cols, value_name
            ))
        }
    } else {
        modi_y <- y
    }
    col_sums <- Matrix::colSums(modi_y)
    col_sums[col_sums == 0] <- NA
    na_cols <- is.na(col_sums)
    inv_col_sums <- 1 / col_sums
    weights <- modi_y %*% Matrix::Diagonal(x = inv_col_sums)
    # return projection result
    result <- as.vector(x %*% weights)
    result[na_cols] <- NA
    result
}
