#' Generate Pseudobulk V(D)J Feature Space
#'
#' This function creates a pseudobulk V(D)J feature space from single-cell data,
#' aggregating V(D)J information into pseudobulk groups. It supports input as
#' either a `Milo` object or a `SingleCellExperiment` object.
#'
#' @param milo A `Milo` or `SingleCellExperiment` object containing V(D)J data.
#' @param pbs Optional. A binary matrix with cells as rows and pseudobulk groups
#'  as columns.
#'  - If `milo` is a `Milo` object, this parameter is not required.
#'  - If `milo` is a `SingleCellExperiment` object, either `pbs` or
#' `col_to_bulk` must be provided.
#' @param col_to_bulk Optional character or character vector. Specifies
#' `colData` column(s) to generate `pbs`. If multiple columns are provided,
#' they will be combined. Default is `NULL`.
#'  - If `milo` is a `Milo` object, this parameter is not required.
#'  - If `milo` is a `SingleCellExperiment` object, either `pbs` or
#' `col_to_bulk` must be provided.
#' @param col_to_take Optional character or vector of characters. Specifies
#' names of colData of milo that need to identify the most common value for each
#'  pseudobulk Default is `NULL`.
#' @param normalise Logical. If `TRUE`, scales the counts of each V(D)J gene
#' group to 1 for each pseudobulk. Default is `TRUE`.
#' @param renormalize Logical. If `TRUE`, rescales the counts of each V(D)J gene
#' group to 1 for each pseudobulk after removing 'missing' calls. Useful when
#' `setupVdjPseudobulk()` was run with `remove_missing = FALSE`.
#' Default is `FALSE`.
#' @param min_count Integer. Sets pseudobulk counts in V(D)J gene groups with
#' fewer than this many non-missing calls to 0. Default is `1`.
#' @param mode_option Character. Specifies the mode for extracting V(D)J genes.
#' Must be one of `c('B', 'abT', 'gdT')`. Default is `'abT'`.
#'   - Note: This parameter is considered only when `extract_cols = NULL`.
#'   - If `NULL`, uses column names such as `v_call_VDJ` instead of
#' `v_call_abT_VDJ`.
#' @param extract_cols Character vector. Specifies column names where V(D)J
#' information is stored. Default is `c('v_call_abT_VDJ_main',
#' 'j_call_abT_VDJ_main', ' 'v_call_abT_VJ_main', 'j_call_abT_VJ_main')`.
#' @param verbose Logical. If `TRUE`, prints messages and warnings.
#' Default is `TRUE`.
#' @details
#' This function aggregates V(D)J data into pseudobulk groups
#' based on the following logic:
#' - **Input Requirements**:
#'  - If `milo` is a `Milo` object, neither `pbs` nor `col_to_bulk` is required.
#'  - If `milo` is a `SingleCellExperiment` object, the user must provide either
#'  `pbs` or `col_to_bulk`.
#' - **Normalization**:
#'  - When `normalise = TRUE`, scales V(D)J counts to 1 for each pseudobulk
#' group.
#'  - When `renormalize = TRUE`, rescales the counts after removing
#' 'missing' calls.
#' - **Mode Selection**:
#'  - If `extract_cols = NULL`, the function relies on `mode_option` to
#' determine which V(D)J columns to extract.
#' - **Filtering**:
#'  - Uses `min_count` to filter pseudobulks with insufficient counts for
#' V(D)J groups.
#'
#' @examples
#' data(sce_vdj)
#' sce_vdj <- setupVdjPseudobulk(sce_vdj,
#'     already.productive = FALSE,
#'     allowed_chain_status = c("Single pair", "Extra pair")
#' )
#' # Build Milo Object
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
#' # Construct pseudobulked VDJ feature space
#' pb.milo <- vdjPseudobulk(milo_object, col_to_take = "anno_lvl_2_final_clean")
#'
#' @return SingleCellExperiment object
#' @include check.R
#' @include getPbs.R
#' @import SingleCellExperiment
#' @importFrom rlang abort
#' @importFrom methods is
#' @importFrom Matrix t
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom stats model.matrix contrasts
#' @importFrom miloR Milo
#'
#' @export
vdjPseudobulk <- function(
    milo, pbs = NULL, col_to_bulk = NULL,
    extract_cols = c(
        "v_call_abT_VDJ_main",
        "j_call_abT_VDJ_main",
        "v_call_abT_VJ_main",
        "j_call_abT_VJ_main"
    ),
    mode_option = c(
        "abT",
        "gdT", "B"
    ),
    col_to_take = NULL,
    normalise = TRUE,
    renormalize = FALSE,
    min_count = 1L,
    verbose = TRUE) {
    if (!is(milo, "Milo") && !is(milo, "SingleCellExperiment")) {
        abort("Must be either Milo or SingleCellExperiment object") # nocov
    }
    .classCheck(pbs, "Matrix")
    if (!all(col_to_bulk %in% names(colData(milo)))) {
        abort("col_to_bulk should within the name of coldata of milo") # nocov
    }
    if (!all(col_to_take %in% names(colData(milo)))) {
        abort("col_to_take should within the name of coldata of milo") # nocov
    }
    .typeCheck(normalise, "logical")
    .typeCheck(renormalize, "logical")
    .typeCheck(min_count, "numeric")
    min_count <- as.integer(min_count)
    .typeCheck(mode_option, "character")
    mode_option <- match.arg(mode_option)
    .typeCheck(extract_cols, "character")
    if (is(milo, "Milo")) {
        pbs <- miloR::nhoods(milo)
    } else {
        pbs <- .getPbs(pbs, col_to_bulk, milo, verbose)
    }
    extract_cols <- .determExtractColN(extract_cols, mode_option, milo)
    pseudo_vdj_feature <- .featureSpaceConstruct(milo, extract_cols, pbs)
    if (normalise) {
        pseudo_vdj_feature <- .normalizeFeatureSpace(
            pseudo_vdj_feature,
            extract_cols, min_count, renormalize, milo
        )
    }
    pb.milo <- .packFeatureSpace(pbs, col_to_take, milo, pseudo_vdj_feature)
    return(pb.milo)
}


#' determine the columns in the colData where the main VDJ information is stored
#'
#' @param extract_cols names of columns in the colData where the main VDJ
#' information stores, passed from vdjPseudobulk
#' @param mode_option  Specifies the mode for extracting V(D)J genes
#' @param milo Milo or SingleCellExperiment object provided by user
#' @keywords internal
#' @return a character vector stores the names of columns in the colData where
#' the main VDJ information stores
#' @import SingleCellExperiment
#' @importFrom S4Vectors DataFrame
.determExtractColN <- function(extract_cols, mode_option, milo) {
    if (is.null(extract_cols)) {
        if (is.null(mode_option)) { # nocov start
            all.col.n <- colnames(colData(milo))
            extract_cols <- all.col.n[grep(
                "_call_VDJ_main|_call_VJ_main",
                all.col.n
            )]
        } else {
            all.col.n <- colnames(colData(milo))
            extract_cols <- all.col.n[grep(paste0(
                mode_option, "_VDJ_main|", mode_option,
                "_call_VJ_main"
            ), all.col.n)]
        } # nocov end
    }
    return(extract_cols)
}


#' Construct VDJ feature space
#'
#' @param milo Milo or SingleCellExperiment object provided by user
#' @param extract_cols columns of names where to extract the VDJ information
#' @keywords internal
#' @import SingleCellExperiment
#' @param pbs cell x pseudobulk adjacent matrix
#' @importFrom S4Vectors DataFrame
#' @importFrom stats model.matrix contrasts
#' @importFrom Matrix t
#' @return constructed feature space
.featureSpaceConstruct <- function(milo, extract_cols, pbs) {
    vjs0 <- data.frame(colData(milo)[extract_cols])
    ## model.matrix need factor input
    vjs0[] <- lapply(vjs0, function(x) {
        if (!is.factor(x)) {
            as.factor(x)
        } else {
            x # nocov
        }
    })
    one_hot_encoded <- model.matrix(~ . - 1,
        data = vjs0,
        contrasts.arg = lapply(vjs0,
            contrasts,
            contrasts = FALSE
        )
    ) # prevent reference level
    colnames(one_hot_encoded) <- gsub(
        "^[^.]*\\main", "",
        colnames(one_hot_encoded)
    )
    pseudo_vdj_feature <- t(t(one_hot_encoded) %*% pbs)
    #  an dgeMatrix with dim pseudobulk x vdj
    return(pseudo_vdj_feature)
}


#' Normalize Feature Space
#'
#' Make sure the sum of each V, D, and J gene within a pseudobulk equals to 1
#'
#' @param pseudo_vdj_feature constructed feature space
#' @param extract_cols names of columns to extract the VDJ information
#' @param min_count the minim count of a V/D/J gene
#' @param renormalize Whether to renormalize the matrix
#' @param milo Milo or SingleCellExperiment object provided by user
#' @keywords internal
#' @return normalized VDJ feature space
.normalizeFeatureSpace <- function(
    pseudo_vdj_feature, extract_cols, min_count,
    renormalize, milo) {
    ## identify any missing calls inserted by the setup, will end with
    ## _missing negate as we want to actually remove them later
    define.mask <- rep(TRUE, length(colnames(pseudo_vdj_feature)))
    define.mask[grep("_missing", colnames(pseudo_vdj_feature))] <- FALSE
    # loop over V(D)J categories
    pseudo_vdj_feature <- Reduce(
        function(psef, col_n) {
            .normalizePerVDJ(
                psef, col_n,
                renormalize = renormalize, define.mask = define.mask,
                milo = milo, min_count = min_count
            )
        }, extract_cols,
        init = pseudo_vdj_feature
    )
    return(pseudo_vdj_feature)
}


#' function to normalize a specific kind of VDJ gene in feature space
#'
#' @param pseudo_vdj_feature constructed feature space
#' @param col_n name of column to extract the VDJ information
#' @param renormalize Whether to renormalize the matrix
#' @param define.mask logical vector determine whether the V/D/J gene
#' should be masked when normalizing
#' @param milo Milo or SingleCellExperiment object provided by user
#' @param min_count the minim count of a V/D/J gene
#' @keywords internal
#' @import SingleCellExperiment
#' @return feature space normalized on specifed V/D/J gene
.normalizePerVDJ <- function(
    pseudo_vdj_feature, col_n, renormalize,
    define.mask, milo, min_count) {
    # identify columns holding genes belonging to the category and then
    # normalise the values to 1 for each pseudobulk
    group.mask <- colnames(pseudo_vdj_feature) %in%
        unique(colData(milo)[[col_n]])
    # identify the defined (non-missing) calls for the group
    group.define.mask <- define.mask & group.mask
    # compute sum of of non-missing values for each pseudobulk for this
    # category and compare to the min_count
    define.count <- apply(pseudo_vdj_feature[, group.define.mask], 1, sum)
    defined.min.counts <- define.count >= min_count
    # normalise the pseudobulks
    pseudo_vdj_feature[, group.define.mask] <- pseudo_vdj_feature[
        ,
        group.define.mask
    ] / define.count
    if (renormalize) { # nocov start
        redefine.count <- apply(
            pseudo_vdj_feature[defined.min.counts, group.define.mask],
            1, sum
        )
        pseudo_vdj_feature[
            defined.min.counts,
            group.define.mask
        ] <- pseudo_vdj_feature[
            defined.min.counts,
            group.define.mask
        ] / define.count
    } # nocov end
    pseudo_vdj_feature[!defined.min.counts, group.define.mask] <- 0
    return(pseudo_vdj_feature)
}

#' Pack the normalized feature space into new Milo object
#'
#' The metadata will derived from the original milo
#' @param pbs cell x pseudobulk adjacent matrix
#' @param col_to_take Optional character or vector of characters.
#' Specifies names of colData of milo that need to identify the most
#' common value for each pseudobulk
#' @param milo Milo or SingleCellExperiment object provided by user
#' @param pseudo_vdj_feature VDJ feature space
#' @keywords internal
#' @return Milo object with VDJ feature space stored in its assay
#' @include getPbs.R
#' @import SingleCellExperiment
#' @importFrom miloR Milo
#' @importFrom rlang abort
#' @importFrom Matrix t
#' @importFrom S4Vectors SimpleList DataFrame
.packFeatureSpace <- function(pbs, col_to_take, milo, pseudo_vdj_feature) {
    # create colData for the new pesudobulk object milo@metadata$feature.space
    # <- pseudo_vdj_feature
    pbs.col <- .getPbsCol(pbs, col_to_take = col_to_take, milo = milo)

    # create a new SingelCellExperiment object as result
    pb.sce <- SingleCellExperiment(
        assay = SimpleList(Feature_space = t(pseudo_vdj_feature)),
        rowData = DataFrame(
            row.names =
                colnames(pseudo_vdj_feature)
        ), colData = pbs.col
    )
    # store pseudobulk assignment, as a sparse for storage efficiency transpose
    # as the original matrix is cells x pseudobulks
    pb.milo <- Milo(pb.sce)
    miloR::nhoods(pb.milo) <- t(pbs)
    return(pb.milo)
}
