#' vdj_pseudobulk
#'
#' making pseudobulk vdj feature space
#' @param milo milo object or SingleCellExperiment object
#' @param pbs matrix
#' Optional binary matrix with cells as rows and pseudobulk groups as columns,
#' if milo is a milo object, no need to provide
#' if milo is a SingleCellExperiment object, user should only provide either pbs or obs_to_bulk
#' @param obs_to_bulk str or a list of str, NULL by default
#' Optional obs column(s) to group pseudobulks into; if multiple are provided, they will be combined
#' if milo is a milo object, no need to provide
#' if milo is a SingleCellExperiment object, user should only provide either pbs or obs_to_bulk
#' @param obs_to_take str or a list of str, NULL by default
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
#' @param extract_cols list of str
#' c('v_call_abT_VDJ_main', 'j_call_abT_VDJ_main', 'v_call_abT_VJ_main', 'j_call_abT_VJ_main') by default
#'  Column names where VDJ/VJ information is stored so that this will be used instead of the standard columns.
#' @return SingleCellEXperiment
#' @export
vdj.pseudobulk <- function(milo, pbs = NULL, obs_to_bulk = NULL, obs_to_take = NULL,
  normalise = TRUE, renormalise = FALSE, min_count = 1L, extract_cols = c("v_call_abT_VDJ_main",
    "j_call_abT_VDJ_main", "v_call_abT_VJ_main", "j_call_abT_VJ_main"), mode_option = c("abT",
    "gdT", "B")) {
  
  # type check
  if (!is(milo, "Milo") && !is(milo, "SingleCellExperiment")) {
    abort("uncompatible data type, \nmilo msut be either Milo or SingleCellExperiment object")
  }
  .class.check(pbs, "Matrix")
  if (!all(obs_to_bulk %in% names(colData(milo_object))))
    abort("Inappropriate argument value: \nobs_to_bulk should within the name of coldata of milo")
  if (!all(obs_to_take %in% names(colData(milo_object))))
    abort("Inappropriate argument value: \nobs_to_take should within the name of coldata of milo")
  .type.check(normalise, "logical")
  .type.check(renormalise, "logical")
  min_count <- as.integer(min_count)
  .type.check(mode_option, "character")
  mode_option <- match.arg(mode_option)
  .type.check(extract_cols, "character")



  # determ the value of pbs
  if (is(milo, "Milo"))
    pbs = milo@nhoods else pbs = .get.pbs(pbs, obs_to_bulk, milo)

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
  one_hot_encoded <- model.matrix(~. - 1, data = vjs0, contrasts.arg = lapply(vjs0,
    contrasts, contrasts = FALSE))  # prevent reference level
  colnames(one_hot_encoded) <- gsub("^[^.]*\\main", "", colnames(one_hot_encoded))
  pseudo_vdj_feature <- t(t(one_hot_encoded) %*% pbs)  #  pseudobulk x vdj
  # this is an dgeMatrix, which cannot coerce to data.frame, go ahead first.
  if (normalise) {
    ## identify any missing calls inserted by the setup, will end with _missing
    ## negate as we want to actually remove them later
    define.mask <- rep(TRUE, length(colnames(pseudo_vdj_feature)))
    define.mask[grep("_missing", colnames(pseudo_vdj_feature))] <- FALSE
    # loop over V(D)J categories
    for (col_n in extract_cols) {
      # identify columns holding genes belonging to the category and then
      # normalise the values to 1 for each pseudobulk
      group.mask <- colnames(pseudo_vdj_feature) %in% unique(milo@colData[[col_n]])

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
  # create obs for the new pesudobulk object milo@metadata$feature.space <-
  # pseudo_vdj_feature
  pbs.obs <- .get.pbs.obs(pbs, obs_to_take = obs_to_take, milo = milo)
  # create a new SingelCellExperiment object as result
  pb.sce <- SingleCellExperiment(assay = SimpleList(X = t(pseudo_vdj_feature)),
    rowData = DataFrame(row.names = colnames(pseudo_vdj_feature)), colData = pbs.obs)
  # store pseudobulk assignment, as a sparse for storage efficiency transpose
  # as the original matrix is cells x pseudobulks

  pb.milo <- Milo(pb.sce)
  nhoods(pb.milo) <- t(nhoods(milo))
  return(pb.milo)
}