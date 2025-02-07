#' .getPbs
#'
#' Helper function to ensure we have cells by pseudobulks matrix which we can use for pseudobulking.
#' @param pbs pbs parameter provided by vdjPseudobulk(),  cells by pseudobulks matrix or NULL
#' @param col_to_bulk col_to_bulk parameter provided by vdjPseudobulk(), column's name of colData from milo
#' @param milo SingleCellExperiment object
#' @param verbose logical, whether to print messages
#' @importFrom rlang abort
#' @importFrom stats model.matrix
#' @importFrom Matrix Matrix
#' @keywords internal
#' @return a cell x pseudobulk matrix
.getPbs <- function(pbs, col_to_bulk, milo, verbose = TRUE) {
    # some way to pseudobulk
    if (is.null(pbs) && is.null(col_to_bulk)) {
        abort("You must specify 'pbs' or 'col_to_bulk'.")
    }
    # but just one
    if (!is.null(pbs) && !is.null(col_to_bulk)) {
        abort("You must specify 'pbs' or 'col_to_bulk', not both.")
    }
    if (!is.null(pbs)) {
        return(pbs)
    }
    if (!is.null(col_to_bulk)) {
        msg <- paste0(col_to_bulk, collapse = ", ")
        if (verbose) {
            message(sprintf("Generating pseudobulks according to colData %s ...", msg),
                appendLF = FALSE
            )
        }
        tobulk <- lapply(col_to_bulk, function(x) {
            colData(milo)[[x]]
        })
        names(tobulk) <- col_to_bulk
        tobulk <- as.data.frame(tobulk)
        tobulk <- as.data.frame(apply(tobulk, 1, paste, collapse = ",", simplify = FALSE))
        tobulk <- model.matrix(~ t(tobulk) - 1)
        colnames(tobulk) <- gsub("t\\(tobulk\\)", "", colnames(tobulk))
        tobulk <- Matrix(tobulk, sparse = TRUE)
        if (verbose) {
            message("Complete")
            message(sprintf("The number of pseudobulks is %d", dim(tobulk)[2]))
        }
        return(tobulk)
    }
}

#' .getPbsCol
#'
#' Helper function to create the new pseudobulk object's coldata.
#' @param pbs dgeMatrix, cell x pseudobulk binary matrix
#' @param col_to_take character vector, names of colData of milo that need to be processed
#' @param milo Milo or SingleCellExperiment object
#' @importFrom Matrix t
#' @importFrom S4Vectors DataFrame
#' @importFrom stats model.matrix contrasts
#' @keywords internal
#' @return pbs_col, a DataFrame which will be passed to the new SingleCellExperiment object as colData of vdj x pseudobulk assays
.getPbsCol <- function(pbs, col_to_take, milo) {
    # prepare per-pseudobulk calls of specified metadata columns
    if (!is.null(col_to_take)) {
        pbs.col <- DataFrame()
        pbs.col <- Reduce(function(pbs.col, anno_col) {
            pbs.col <- .getPbsPerCol(pbs.col, anno_col, milo, pbs)
            return(pbs.col)
        }, col_to_take, pbs.col)
    }
    # report the number of cells for each pseudobulk ensure pbs is an array so
    # that it sums into a vector that can go in easily
    pbs.col["cell_count"] <- apply(pbs, 2, sum)
    return(pbs.col)
}



#' .getPbsPerCol
#'
#' function used in Reduce to get the PbsCol
#' @param pbs.col DataFrame object used to store the result of each iteration
#' @param anno_col colname where to generate the metadata from
#' @param milo milo or SingleCellExperiment objects provided by user
#' @param pbs dgeMatrix, cell x pseudobulk binary matrix
#' @importFrom Matrix t
#' @importFrom S4Vectors DataFrame
#' @importFrom stats model.matrix contrasts
#' @keywords internal
#' @return DataFrame object, serve as part of metadata of the new milo object
.getPbsPerCol <- function(pbs.col, anno_col, milo, pbs) {
    fa <- as.factor(colData(milo)[[anno_col]])
    fa <- data.frame(fa)
    anno.dummies <- model.matrix(~ . - 1,
        data = fa,
        contrasts.arg = lapply(fa,
            contrasts,
            contrasts = FALSE
        )
    )
    anno.count <- Matrix::t(pbs) %*% anno.dummies
    # t(cell x  pseudo)   %*% (cell x element) = pseudo x element
    anno.sum <- apply(anno.count, 1, sum)
    anno.frac <- anno.count / anno.sum
    anno.frac <- as.matrix(anno.frac)
    anno.frac <- DataFrame(anno.frac)
    colnames(anno.frac) <- colnames(anno.dummies)
    pbs.col[anno_col] <- colnames(anno.frac)[apply(anno.frac, 1, which.max)]
    pbs.col[paste0(anno_col, "_fraction")] <- apply(anno.frac, 1, max)
    return(pbs.col)
}
