#' minmax.scale
#'
#' scale the value to range 0 to 1
#' @param data dataframe need to be scale
#' @returns scaled value
.minmax.scale <- function(data) {
    maxv <- apply(data, 2, max)
    minv <- apply(data, 2, min)
    data_scale <- apply(data, 1, `-`, minv)/(maxv - minv)
    t(data_scale)
}
