#' .max.min.sampling
#'
#' function for max min sampling of waypoints
#'
#' @param data data matrix along which to sample the waypoints, usually diffusion components
#' @param num_waypoints number of waypoints to sample
#' @param seed Random number generator seed to find initial guess.
#' @return Series reprenting the sampled waypoints
.max.min.sampling <- function(data, num_waypoints, seed = NULL) {
  message("Sampling and flocking waypoints...")
  if (!is.null(seed))
    set.seed(seed)
  no.iterations <- as.integer(num_waypoints/ncol(data))
  waypoints <- c()
  # Sample along each componet
  N = nrow(data)
  for (ind in colnames(data)) {
    vecs <- data[, ind]
    iter.set <- sample(1:N, 1)
    dists <- matrix(0, nrow = N, ncol = no.iterations)
    dists[, 1] <- abs(vecs - vecs[iter.set])
    for (k in 2:no.iterations - 1) {
      # minium distances
      min_dists <- apply(dists[, 1:k, drop = FALSE], 1, min)
      # new waypoint: where the maxium of the minimum distane locate
      new.wp <- which.max(min_dists)
      iter.set <- c(iter.set, new.wp)
      dists[, k + 1] <- abs(vecs - data[new.wp, ind])
    }
    waypoints <- c(waypoints, iter.set)
  }
  waypoints <- unique(waypoints)
}