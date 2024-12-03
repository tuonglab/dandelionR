#' .get.pbs
#'
#' Helper function to ensure we have cells by pseudobulks matrix which we can use for pseudobulking.
#' @param pbs pbs parameter provided by vdj_pseudobulk(),  cells by pseudobulks matrix or NULL
#' @param col_to_bulk col_to_bulk parameter provided by vdj_pseudobulk(), column's name of colData from milo
#' @param milo SingleCellExperiment object
#' @return a cell x pseudobulk matrix
.get.pbs <- function(pbs, col_to_bulk, milo) {
    # some way to pseudobulk
    requireNamespace("rlang")
    if (is.null(pbs) && is.null(col_to_bulk)) {
        rlang::abort("You must specify 'pbs' or 'col_to_bulk'.")
    }
    # but just one
    if (!is.null(pbs) && !is.null(col_to_bulk)) {
        rlang::abort("You must specify 'pbs' or 'col_to_bulk', not both.")
    }
    if (!is.null(pbs)) {
        return(pbs)
    }
    if (!is.null(col_to_bulk)) {
        msg <- paste0(col_to_bulk, collapse = ", ")
        message(sprintf("Generating pseudobulks according to colData %s ...", msg),
            appendLF = FALSE)
        tobulk <- lapply(col_to_bulk, function(x) {
            colData(milo)[[x]]
        })
        names(tobulk) <- col_to_bulk
        tobulk <- as.data.frame(tobulk)
        tobulk <- as.data.frame(apply(tobulk, 1, paste, collapse = ",", simplify = FALSE))
        requireNamespace("stats")
        tobulk <- stats::model.matrix(~t(tobulk) - 1)
        colnames(tobulk) <- gsub("t\\(tobulk\\)", "", colnames(tobulk))
        requireNamespace("Matrix")
        tobulk <- Matrix::Matrix(tobulk, sparse = TRUE)
        message("Complete")
        message(sprintf("The number of pseudobulks is %d", dim(tobulk)[2]))
        return(tobulk)
    }
}

#' .get.pbs.col
#'
#' Helper function to create the new pseudobulk object's coldata.
#' @param pbs dgeMatrix, cell x pseudobulk binary matrix
#' @param col_to_take character vector, names of colData of milo that need to be processed  
#' @param milo Milo or SingleCellExperiment object
#' @return pbs_col, a DataFrame which will be passed to the new SingleCellExperiment object as colData of vdj x pseudobulk assays
.get.pbs.col <- function(pbs, col_to_take, milo) {
    # prepare per-pseudobulk calls of specified metadata columns
    if (!is.null(col_to_take)) {
        requireNamespace("S4Vectors")
        pbs.col <- S4Vectors::DataFrame()
        for (anno_col in col_to_take) {
            fa <- as.factor(colData(milo)[[anno_col]])
            fa <- data.frame(fa)
            requireNamespace("stats")
            anno.dummies <- stats::model.matrix(~. - 1, data = fa, contrasts.arg = lapply(fa,
                stats::contrasts, contrasts = FALSE))
            requireNamespace("Matrix")
            anno.count <- Matrix::t(pbs) %*% anno.dummies  # t(cell x  pseudo)   %*% (cell x element) = pseudo x element
            anno.sum <- apply(anno.count, 1, sum)
            anno.frac <- anno.count/anno.sum
            anno.frac <- S4Vectors::DataFrame(anno.frac)
            colnames(anno.frac) <- colnames(anno.dummies)
            pbs.col[anno_col] <- colnames(anno.frac)[apply(anno.frac, 1, which.max)]
            pbs.col[paste0(anno_col, "_fraction")] <- apply(anno.frac, 1, max)
        }
    }
    # report the number of cells for each pseudobulk ensure pbs is an array so
    # that it sums into a vector that can go in easily
    pbs.col["cell_count"] <- apply(pbs, 2, sum)
    return(pbs.col)
}
