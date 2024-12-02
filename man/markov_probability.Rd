% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/markov_probability.R
\name{markov_probability}
\alias{markov_probability}
\title{markov_probability}
\usage{
markov_probability(
  milo,
  diffusionmap,
  terminal_state,
  root_cell,
  pseudotime_key,
  scale_components = TRUE,
  num_waypoints = 500
)
}
\arguments{
\item{milo}{milo or SingelCellExperiment object, used to store the result}

\item{diffusionmap}{DiffusionMap object corresponds to milo}

\item{terminal_state}{the index of the terminal state}

\item{root_cell}{the index of the root state}

\item{pseudotime_key}{the column name in the colData that holds the inferred pseudotime}

\item{scale_components}{logical, If True, the components will be scale before constructing markov chain}

\item{num_waypoints}{integer, 500L by default. Number of waypoints to sample to construct markov chain.}
}
\value{
milo or SinglCellExperiment object with pseudotime, probabilities in its colData
}
\description{
Preprocessing data and Construct markov chain and calculate probabilities
}
