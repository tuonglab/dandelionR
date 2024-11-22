#' .class.check
#'
#' check whether the input is with the correct class
#' @param input the input need to be check
#' @param must the type we need
.class.check <- function(input, must) {
  requireNamespace("rlang")
  requireNamespace("methods")
  if (is.null(input))
    return()
  if (!methods::is(input, must)) {
    rlang::abort(paste0("The '", as.character(substitute(input)), "' must be ", must,
      ", not ", class(input)))
  }
}


#' .type.check
#'
#' check whether the input has the correct type
#' @param input the input need to be check
#' @param must the type we need
.type.check <- function(input, must) {
  requireNamespace("methods")
  if (is.null(input))
    return()
  if (!methods::is(input, must)) {
    requireNamespace("rlang")
    requireNamespace("BiocGenerics")
    rlang::abort(paste0("The '", as.character(substitute(input)), "' must be ", must,
      ", not ", BiocGenerics::type(input)))
  }
}