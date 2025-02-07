#' .maxMinSampling
#'
#' function for max min sampling of waypoints
#'
#' @param data data matrix along which to sample the waypoints, usually diffusion components
#' @param num_waypoints number of waypoints to sample
#' @param verbose logical, print progress
#' @keywords internal
#' @return Series reprenting the sampled waypoints
.maxMinSampling <- function(data, num_waypoints, verbose = TRUE) {
    if (verbose) message("Sampling and flocking waypoints...")

    no.iterations <- as.integer(num_waypoints / ncol(data))
    waypoints <- c()
    # Sample along each componet
    waypoints <- Reduce(function(waypoints, ind){
      waypoints <- .waypiontsPerCol(waypoints = waypoints, ind = ind, data = data)
      return(waypoints)
    }, colnames(data), waypoints)
    waypoints <- unique(waypoints)
    return(waypoints)
}


#' find the waypoints according to certain columns of data
#' 
#' @param waypoints integer vector used to store waypoints
#' @param ind columns' colnames
#' @param data scaled diffusionmap
#' @keywords internal
#' @return a numeric vector containing waypoints' index
.waypiontsPerCol <- function(waypoints, ind, data)
{
  N <- nrow(data)
  vecs <- data[, ind]
  iter.set <- sample(seq_len(N), 1)
  dists <- matrix(0, nrow = N, ncol = no.iterations)
  dists[, 1] <- abs(vecs - vecs[iter.set])
  iterdists <- list(iter.set = iter.set, dists = dists)
  iterdists<- Reduce(function(iterdists,k){
    .findNewWaypoints(iterdists, k, vec = vecs, ind = ind)
  }, seq_len(no.iterations-1), iterdists)
  waypoints <- c(waypoints, iterdists$iter.set)
  return(waypoints)
}


#' function used in Reduce to find new waypoint in an iteration
#' 
#' @param iterdists a list containing both waypoints deteted in
#'  the former iterations and the distance matrix used to find waypoints
#' @param k the iteration number
#' @param vec a numeric vector used to calculate distance of waypoints to
#'  each points
#' @param ind colnames
#' @keywords internal
#' @return a list containing updated distance matrix and new waypoints
.findNewWaypoints <- function(iterdists, k, vecs, ind)
{
  iter.set <- iterdists$iter.set
  dists <- iterdists$dists
  # minium distances
  min_dists <- apply(dists[, seq_len(k), drop = FALSE], 1, min)
  # new waypoint: where the maxium of the minimum distane locate
  new.wp <- which.max(min_dists)
  iter.set <- c(iter.set, new.wp)
  dists[, k + 1] <- abs(vecs - data[new.wp, ind])
  return(list(iter.set = iter.set, dists = dists))
}

