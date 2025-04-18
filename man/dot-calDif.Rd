% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/projectProbability.R
\name{.calDif}
\alias{.calDif}
\title{function help to calculate the diffusion distance}
\usage{
.calDif(idx, eigenvector, lambda_t, K)
}
\arguments{
\item{idx}{integer the index of the calculated value}

\item{eigenvector}{numeric vector, the eigenvector from diffusion map}

\item{lambda_t}{eigenvalues to the power of t(diffusion time)}

\item{K}{The number of the eigenvectors to be used in calculation}
}
\value{
updated diffusion distance matrix after one iteration
}
\description{
function help to calculate the diffusion distance
}
\keyword{internal}
