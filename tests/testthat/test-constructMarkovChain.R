library(igraph)
library(Matrix)
library(bluster)
library(purrr)
library(stats)

# Mock data for testing
set.seed(100)
wp_data <- matrix(rnorm(100), ncol = 2)
knn <- 5
pseudotime <- runif(50)
waypoints <- sample(1:50, 10)

# Test for .constructMarkovChain
test_that(".constructMarkovChain works correctly", {
    result <- .constructMarkovChain(wp_data, knn, pseudotime, waypoints)

    expect_true(is(result, "dgCMatrix"))
    expect_equal(nrow(result), length(waypoints))
    expect_equal(ncol(result), length(waypoints))
    expect_true(all(result >= 0))
    expect_true(all(rowSums(result) <= 1))
})

# Additional tests for edge cases
test_that(".constructMarkovChain handles edge cases", {
    # Test with minimum waypoints
    waypoints_min <- sample(1:50, 2)
    result_min <- .constructMarkovChain(wp_data, knn, pseudotime, waypoints_min)
    expect_true(is(result_min, "dgCMatrix"))
    expect_equal(nrow(result_min), length(waypoints_min))
    expect_equal(ncol(result_min), length(waypoints_min))

    # Test with maximum waypoints
    waypoints_max <- sample(1:50, 50)
    result_max <- .constructMarkovChain(wp_data, knn, pseudotime, waypoints_max)
    expect_true(is(result_max, "dgCMatrix"))
    expect_equal(nrow(result_max), length(waypoints_max))
    expect_equal(ncol(result_max), length(waypoints_max))
})
