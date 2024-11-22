#' vdj_pseudobulk
#'
#' making pseudobulk vdj feature space
#' @param milo milo object or SingleCellExperiment object
#' @param pbs Optional binary matrix with cells as rows and pseudobulk groups as columns,
#'  - if milo is a milo object, no need to provide
#'  - if milo is a SingleCellExperiment object, user should only provide either pbs or col_to_bulk
#' @param col_to_bulk character or character vector, NULL by default
#' Optional colData column(s) to generate pbs, if multiple are provided, they will be combined
#'  - if milo is a milo object, no need to provide
#'  - if milo is a SingleCellExperiment object, user should only provide either pbs or col_to_bulk
#' @param col_to_take str or a list of str, NULL by default
#' Optional obs column(s) to identify the most common value of for each pseudobulk.
#' @param normalise bool, True by default
#'  If True, will scale the counts of each V(D)J gene group to 1 for each pseudobulk.
#' @param renormalise bool, False by default
#' If True, will re-scale the counts of each V(D)J gene group to 1 for each pseudobulk with any 'missing' calls removed.
#' Relevant with normalise as True, if setup_vdj_pseudobulk() was ran with remove_missing set to False.
#' @param min_count int, 1 by default
#' Pseudobulks with fewer than these many non-'missing' calls in a V(D)J gene group will have their non-'missing' calls set to 0 for that group. Relevant with normalise as True.
#' @param mode_option, must be one element of the vector c('B','abT','gdT'), 'abT' by default
#' Note: only when you set extract_cols to NULL, will this argument be considered!
#' Optional mode for extracting the V(D)J genes. If set as NULL, it will use e.g. v_call_VD` instead of v_call_abT_VDJ.
#' @param extract_cols character vector
#' with default value c('v_call_abT_VDJ_main', 'j_call_abT_VDJ_main', 'v_call_abT_VJ_main', 'j_call_abT_VJ_main')
#'  Column names where VDJ/VJ information is stored so that this will be used instead of the standard columns.
#' @return SingleCellExperiment object ...
#' @include check.R
#' @include get.pbs.R
#' @import SingleCellExperiment
#' @import miloR
#' @export
vdj_pseudobulk <- function(
    milo, 
    pbs = NULL, 
    col_to_bulk = NULL, 
    extract_cols = c("v_call_abT_VDJ_main", "j_call_abT_VDJ_main", "v_call_abT_VJ_main", "j_call_abT_VJ_main"), 
    mode_option = c("abT", "gdT", "B"),
    col_to_take = NULL,
    normalise = TRUE, 
    renormalise = FALSE, 
    min_count = 1L
    ) {
  # type check
  requireNamespace("methods")
  requireNamespace("rlang")
  if (!methods::is(milo, "Milo") && !methods::is(milo, "SingleCellExperiment")) {
    rlang::abort("Uncompatible data type, \nmilo msut be either Milo or SingleCellExperiment object")
  }
  .class.check(pbs, "Matrix")
  if (!all(col_to_bulk %in% names(colData(milo))))
    rlang::abort("Inappropriate argument value: \nocol_to_bulk should within the name of coldata of milo")
  if (!all(col_to_take %in% names(colData(milo))))
    rlang::abort("Inappropriate argument value: \ncol_to_take should within the name of coldata of milo")
  .type.check(normalise, "logical")
  .type.check(renormalise, "logical")
  .type.check(min_count, "numeric")
  min_count <- as.integer(min_count)
  .type.check(mode_option, "character")
  mode_option <- match.arg(mode_option)
  .type.check(extract_cols, "character")
  # determ the value of pbs
  if (methods::is(milo, "Milo"))
    pbs <- miloR::nhoods(milo) else pbs <- .get.pbs(pbs, col_to_bulk, milo)
  # set the column used in caculation
  if (is.null(extract_cols)) {
    if (is.null(mode_option)) {
      all.col.n <- colnames(colData(milo))
      extract_cols <- all.col.n[grep("_call_VDJ_main|_call_VJ_main", all.col.n)]
    } else {
      all.col.n <- colnames(colData(milo))
      extract_cols <- all.col.n[grep(paste0(mode_option, "_VDJ_main|", mode_option,
        "_call_VJ_main"), all.col.n)]
    }
  }
  # perform matrix multiplication of pseudobulks by cell matrix by a cells by
  # VJs matrix strat off by creating the cell by VJs matix skip the prefix
  # stuff as the VJ genes will be unique in the columns
  vjs0 = data.frame(colData(milo)[extract_cols])
  ## model.matrix need factor input
  vjs0[] <- lapply(vjs0, function(x) if (!is.factor(x))
    as.factor(x) else x)
  requireNamespace("stats")
  one_hot_encoded <- stats::model.matrix(~. - 1, data = vjs0, contrasts.arg = lapply(vjs0,
    stats::contrasts, contrasts = FALSE))  # prevent reference level
  colnames(one_hot_encoded) <- gsub("^[^.]*\\main", "", colnames(one_hot_encoded))
  requireNamespace("Matrix")
  pseudo_vdj_feature <- Matrix::t(t(one_hot_encoded) %*% pbs) #  an dgeMatrix with dim pseudobulk x vdj
  if (normalise) {
    ## identify any missing calls inserted by the setup, will end with _missing
    ## negate as we want to actually remove them later
    define.mask <- rep(TRUE, length(colnames(pseudo_vdj_feature)))
    define.mask[grep("_missing", colnames(pseudo_vdj_feature))] <- FALSE
    # loop over V(D)J categories
    for (col_n in extract_cols) {
      # identify columns holding genes belonging to the category and then
      # normalise the values to 1 for each pseudobulk
      group.mask <- colnames(pseudo_vdj_feature) %in% unique(colData(milo)[[col_n]])
      # identify the defined (non-missing) calls for the group
      group.define.mask <- define.mask & group.mask
      # compute sum of of non-missing values for each pseudobulk for this
      # category and compare to the min_count
      
      define.count <- apply(pseudo_vdj_feature[, group.define.mask], 1, sum)
      defined.min.counts <- define.count >= min_count
      # normalise the pseudobulks
      pseudo_vdj_feature[, group.define.mask] <- pseudo_vdj_feature[, group.define.mask]/define.count
      if (renormalise) {
        redefine.count <- apply(pseudo_vdj_feature[defined.min.counts, group.define.mask],
          1, sum)
        pseudo_vdj_feature[defined.min.counts, group.define.mask] <- pseudo_vdj_feature[defined.min.counts,
          group.define.mask]/define.count
      }
      pseudo_vdj_feature[!defined.min.counts, group.define.mask] <- 0
    }
  }
  # create colData for the new pesudobulk object milo@metadata$feature.space <-
  # pseudo_vdj_feature
  pbs.col <- .get.pbs.col(pbs, col_to_take = col_to_take, milo = milo)
  # create a new SingelCellExperiment object as result
  requireNamespace("S4Vectors")
  pb.sce <- SingleCellExperiment(assay = S4Vectors::SimpleList(X = Matrix::t(pseudo_vdj_feature)),
    rowData = S4Vectors::DataFrame(row.names = colnames(pseudo_vdj_feature)), colData = pbs.col)
  # store pseudobulk assignment, as a sparse for storage efficiency transpose
  # as the original matrix is cells x pseudobulks
  pb.milo <- Milo(pb.sce)
  nhoods(pb.milo) <- Matrix::t(pbs)
  return(pb.milo)
}