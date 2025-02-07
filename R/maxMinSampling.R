#' .maxMinSampling
#'
#' function for max min sampling of waypoints
#'
#' @param datas data matrix along which to sample the waypoints, usually
#'  diffusion components
#' @param num_waypoints number of waypoints to sample
#' @param verbose logical, print progress
#' @keywords internal
#' @return Series reprenting the sampled waypoints
.maxMinSampling <- function(datas, num_waypoints, verbose) {
    if (verbose) message("Sampling and flocking waypoints...")
    no.iterations <- as.integer(num_waypoints / ncol(datas))
    waypoints <- c()
    # Sample along each componet
    waypoints <- Reduce(function(waypoints, ind) {
        waypoints <- .waypiontsPerCol(
            waypoints = waypoints, ind = ind,
            datas = datas, no.iterations = no.iterations
        )
        return(waypoints)
    }, colnames(datas), waypoints)
    waypoints <- unique(waypoints)
    return(waypoints)
}


#' find the waypoints according to certain columns of data
#'
#' @param waypoints integer vector used to store waypoints
#' @param ind columns' colnames
#' @param datas scaled diffusionmap
#' @keywords internal
#' @return a numeric vector containing waypoints' index
.waypiontsPerCol <- function(waypoints, ind, datas, no.iterations) {
    N <- nrow(datas)
    vecs <- datas[, ind]
    iter.set <- sample(seq_len(N), 1)
    dists <- matrix(0, nrow = N, ncol = no.iterations)
    dists[, 1] <- abs(vecs - vecs[iter.set])
    iterdists <- list(iter.set = iter.set, dists = dists)
    iterdists <- Reduce(function(iterdists, k) {
        .findNewWaypoints(iterdists, k, vecs = vecs, ind = ind, datas = datas)
    }, seq_len(no.iterations - 1), iterdists)
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
.findNewWaypoints <- function(iterdists, k, vecs, ind, datas) {
    iter.set <- iterdists$iter.set
    dists <- iterdists$dists
    # minium distances
    min_dists <- apply(dists[, seq_len(k), drop = FALSE], 1, min)
    # new waypoint: where the maxium of the minimum distane locate
    new.wp <- which.max(min_dists)
    iter.set <- c(iter.set, new.wp)
    dists[, k + 1] <- abs(vecs - datas[new.wp, ind])
    return(list(iter.set = iter.set, dists = dists))
}
