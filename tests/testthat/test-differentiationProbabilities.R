library(testthat)
library(Matrix)
library(MASS)

# Mock data for testing
set.seed(100)
wp_data <- matrix(rnorm(100), ncol = 2)
knn <- 5
pseudotime <- runif(50)
waypoints <- sample(1:50, 10)
terminal_states <- sample(waypoints, 2)

# Test for differentiationProbabilities
test_that("differentiationProbabilities works correctly", {
    result <- differentiationProbabilities(wp_data, terminal_states, knn, pseudotime, waypoints)

    expect_true(is.matrix(result) || is.data.frame(result))
    expect_equal(nrow(result), length(waypoints))
    expect_equal(ncol(result), length(terminal_states))
    expect_true(all(result >= 0))
    expect_true(all(result <= 1))
})

# Test with default terminal_states
test_that("differentiationProbabilities works with default terminal_states", {
    result <- differentiationProbabilities(wp_data, knn = knn, pseudotime = pseudotime, waypoints = waypoints)

    expect_true(is.matrix(result) || is.data.frame(result))
    expect_equal(nrow(result), length(waypoints))
    expect_true(all(result >= 0))
    expect_true(all(result <= 1))
})

# Additional tests for edge cases
test_that("differentiationProbabilities handles edge cases", {
    # Test with minimum waypoints
    waypoints_min <- sample(1:50, 2)
    result_min <- differentiationProbabilities(wp_data, terminal_states, knn, pseudotime, waypoints_min)
    expect_true(is.matrix(result_min) || is.data.frame(result_min))
    expect_equal(nrow(result_min), length(waypoints_min))

    # Test with maximum waypoints
    waypoints_max <- sample(1:50, 50)
    result_max <- differentiationProbabilities(wp_data, terminal_states, knn, pseudotime, waypoints_max)
    expect_true(is.matrix(result_max) || is.data.frame(result_max))
    expect_equal(nrow(result_max), length(waypoints_max))
})
