#' .classCheck
#'
#' check whether the input is with the correct class
#' @param input the input need to be check
#' @param must the type we need
#' @keywords internal
#' @importFrom rlang abort
#' @importFrom methods is
#' @return whether or not the input is the correct class
.classCheck <- function(input, must) {
    if (is.null(input)) {
        return()
    }
    if (!is(input, must)) {
        abort(sprintf(
            "The '%s' must be %s, not %s",
            as.character(substitute(input)), must, class(input)
        ))
    }
}


#' .typeCheck
#'
#' check whether the input has the correct type
#' @param input the input need to be check
#' @param must the type we need
#' @keywords internal
#' @importFrom rlang abort
#' @importFrom methods is
#' @importFrom BiocGenerics type
#' @return whether or not the input is the correct type
.typeCheck <- function(input, must) {
    if (is.null(input)) {
        return()
    }
    if (!is(input, must)) {
        abort(sprintf(
            "The '%s' must be %s, not %s",
            as.character(substitute(input)), must, type(input)
        ))
    }
}
