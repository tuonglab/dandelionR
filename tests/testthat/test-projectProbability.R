library(testthat)
library(destiny)

# Mock data for testing
set.seed(100)
data <- matrix(rnorm(100), ncol = 10)
dm <- DiffusionMap(data)
waypoints <- sample(1:10, 5)
probabilities <- runif(5)

# Test for projectProbability
test_that("projectProbability works correctly", {
    result <- projectProbability(dm, waypoints, probabilities)

    expect_true(is.matrix(result) || is.data.frame(result))
    expect_equal(nrow(result), nrow(data))
    expect_equal(ncol(result), 1)
    expect_true(all(result >= 0))
    expect_true(all(result <= 1))
})

# Additional tests for edge cases
test_that("projectProbability handles edge cases", {
    # Test with minimum waypoints
    waypoints_min <- sample(1:10, 1)
    probabilities_min <- runif(1)
    result_min <- projectProbability(dm, waypoints_min, probabilities_min)
    expect_true(is.matrix(result_min) || is.data.frame(result_min))
    expect_equal(nrow(result_min), nrow(data))
    expect_equal(ncol(result_min), 1)

    # Test with maximum waypoints
    waypoints_max <- sample(1:10, 10)
    probabilities_max <- runif(10)
    result_max <- projectProbability(dm, waypoints_max, probabilities_max)
    expect_true(is.matrix(result_max) || is.data.frame(result_max))
    expect_equal(nrow(result_max), nrow(data))
    expect_equal(ncol(result_max), 1)
})
