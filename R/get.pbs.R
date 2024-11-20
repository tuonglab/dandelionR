#' .get.pbs
#'
#' Helper function to ensure we have cells by pseudobulks matrix which we can use for pseudobulking.
#' @param pbs pbs parameter provided by vdj_pseudobulk(),  cells by pseudobulks matrix or NULL
#' @param col_to_bulk col_to_bulk parameter provided by vdj_pseudobulk(), column's name of colData from milo
#' @param milo SingleCellExperiment object
#'
.get.pbs <- function(pbs, col_to_bulk, milo) {
  # some way to pseudobulk
  if (is.null(pbs) && is.null(obs_to_bulk)) {
    abort("You must specify 'pbs' or 'obs_to_bulk'.")
  }
  # but just one
  if (!is.null(pbs) && !is.null(obs_to_bulk)) {
    abort("You must specify 'pbs' or 'obs_to_bulk', not both.")
  }
  if (!is.null(pbs)) {
    return(pbs)
  }
  if (!is.null(obs_to_bulk)) {
    message(paste0("Generating pseudobulks according to colData ",paste(obs_to_bulk,collapse = ", "),"..."),appendLF = FALSE)
    tobulk <- lapply(obs_to_bulk, function(x) {
      colData(milo)[[x]]
    })
    names(tobulk) <- obs_to_bulk
    tobulk <- as.data.frame(tobulk)
    tobulk <- as.data.frame(apply(tobulk, 1, paste, collapse = ",", simplify = FALSE))
    tobulk <- model.matrix(~t(tobulk) - 1)
    colnames(tobulk) <- gsub("t\\(tobulk\\)", "", colnames(tobulk))
    tobulk <- Matrix(tobulk, sparse = TRUE)
    message("Complete")
    message(paste("The number of pseudobulk is",dim(tobulk)[2]))
    return(tobulk)
  }
}