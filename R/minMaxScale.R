#' minMaxScale
#'
#' scale the value to range 0 to 1
#' @param data dataframe need to be scale
#' @keywords internal
#' @returns scaled value
.minMaxScale <- function(data) {
    maxv <- apply(data, 2, max)
    minv <- apply(data, 2, min)
    data_scale <- apply(data, 1, `-`, minv) / (maxv - minv)
    return(t(data_scale))
}
