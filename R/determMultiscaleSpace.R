#' .determineMultiscaleSpace
#'
#' @param diffusionmap DiffusionMap object
#' @param n_eigs integer, default is NULL. Number of eigen vectors to use.
#' - If is not specified, the number of eigen vectors will be determined using
#'  the eigen gap.
#' @keywords internal
#' @importFrom destiny eigenvectors eigenvalues
#' @returns dataframe
.determineMultiscaleSpace <- function(diffusionmap, n_eigs = NULL) {
    .classCheck(diffusionmap, "DiffusionMap")
    eigenvect <- eigenvectors(diffusionmap)
    eigenval <- eigenvalues(diffusionmap)
    # determine n_eigs
    if (is.null(n_eigs)) {
        val_gaps <- eigenval[2:length(eigenval) - 1] -
            eigenval[2:length(eigenval)]
        n_eigs <- order(val_gaps)[length(val_gaps)]
        if (n_eigs < 3) {
            n_eigs <- order(val_gaps)[length(val_gaps) - 1]
        }
    }
    # Scale the data
    use_eigs <- seq(1, n_eigs)
    eigenvect <- eigenvect[, use_eigs]
    eigenval <- eigenval[use_eigs]
    multis <- apply(eigenvect, 1, `*`, (eigenval / (1 - eigenval)))
    return(data.frame(t(multis)))
}
