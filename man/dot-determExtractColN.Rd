% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vdjPseudobulk.R
\name{.determExtractColN}
\alias{.determExtractColN}
\title{determine the columns in the colData where the main VDJ information is stored}
\usage{
.determExtractColN(extract_cols, mode_option, milo)
}
\arguments{
\item{extract_cols}{names of columns in the colData where the main VDJ
information stores, passed from vdjPseudobulk}

\item{mode_option}{Specifies the mode for extracting V(D)J genes}

\item{milo}{Milo or SingleCellExperiment object provided by user}
}
\value{
a character vector stores the names of columns in the colData where
the main VDJ information stores
}
\description{
determine the columns in the colData where the main VDJ information is stored
}
\keyword{internal}
