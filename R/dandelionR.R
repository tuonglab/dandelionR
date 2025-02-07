#' dandelionR: Single-cell immune repertoire trajectory analysis
#'
#' dandelionR is an R package for performing single-cell immune repertoire
#' trajectory analysis, based on the original python implementation.
#' It provides the necessary functions to interface with scRepertoire and a
#' custom implementation of an absorbing Markov chain for pseudotime
#' inference, inspired by the Palantir Python package.
#'
#' @section Main functions:
#' - \code{\link{setupVdjPseudobulk}}: Preprocess V(D)J Data for Pseudobulk
#' Analysis.
#' - \code{\link{vdjPseudobulk}}: Generate Pseudobulk V(D)J Feature Space.
#' - \code{\link{markovProbability}}: Markov Chain Construction and
#' Probability Calculation.
#' - \code{\link{projectPseudotimeToCell}}: Project Pseudotime and Branch
#' Probabilities to Single Cells.
#'
#' @section Vignettes:
#' See the package vignettes for detailed workflows:
#' \code{vignette('dandelionR')}
#'
#' @section Installation:
#' To install from Bioconductor, use:
#' \preformatted{
#' if (!requireNamespace('BiocManager', quietly = TRUE))
#'     install.packages('BiocManager')
#' BiocManager::install('dandelionR')
#' }
#'
#' @name dandelionR
#' @keywords internal
"_PACKAGE"
