#' .filterCells
#'
#' Helper function that identifies filter_pattern hits in determined column of sce, and then either removes the offeending cells or masks the matched values with a uniform value of '(column's name)_missing'
#' @param sce SingleCellExperiment object, adata in python
#' data after combineTCR, contain both vdj and seq
#' @param col_n mode for extraction the V(D)J genes.
#' @param filter_pattern character string, optional ',|None|No_contig' by default
#' @param remove_missing bool, True by default
#'  - If `TRUE`, will remove cells with contigs matching the filter from the object.
#'  - If `FALSE`, will mask them with a uniform value dependent on the column name.
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData<-
#' @keywords internal
#' @importFrom rlang abort
#' @return filtered SingleCellExperiment object according to the parameter.
.filterCells <- function(sce, col_n, filter_pattern = ",|None|No_contig", remove_missing = TRUE) {
    if (ncol(sce) < 1) {
        abort("None column remains, please check whether the filtering option is correct.")
    }
    # find filter pattern hits in our column of interest
    temp <- colData(sce)[col_n]
    temp <- as.data.frame(temp@listData)
    mask <- grep(filter_pattern, temp[, 1])
    if (remove_missing) {
        sce <- sce[, setdiff(seq_len(ncol(sce)), mask)]
    } else {
        # uniformly mask the filter pattern hits
        colData(sce)[mask, col_n] <- sprintf("%s_missing", col_n)
    }
    return(sce)
}
