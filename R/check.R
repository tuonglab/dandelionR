#' .class.check
#'
#' check whether the input is with the correct class
#' @param input the input need to be check
#' @param must the type we need
#' @import methods
#' @import rlang
.class.check <- function(input, must) {
  if (is.null(input))
    return()
  if (!is(input, must)) {
    abort(paste0("The '", as.character(substitute(input)), "' must be ", must,
      ", not ", class(input)))
  }
}


#' .type.check
#'
#' check whether the input has the correct type
#' @param input the input need to be check
#' @param must the type we need
#' @import methods
#' @import rlang
#' @import BiocGenerics
.type.check <- function(input, must) {
  if (is.null(input))
    return()
  if (!is(input, must)) {
    abort(paste0("The '", as.character(substitute(input)), "' must be ", must,
      ", not ", type(input)))
  }
}