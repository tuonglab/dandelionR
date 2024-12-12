library(igraph)
library(purrr)
library(spam)
library(Matrix)

# Mock data for testing
set.seed(100)
wp_data <- matrix(rnorm(100), ncol = 2)
pseudotime <- runif(50)
waypoints <- sample(seq_len(50), 10)
Transmat <- matrix(runif(100), nrow = 10, ncol = 10)
Transmat <- Transmat / rowSums(Transmat)

# Test for .terminalStateFromMarkovChain
test_that(".terminalStateFromMarkovChain works correctly", {
    result <- .terminalStateFromMarkovChain(Transmat, wp_data, pseudotime, waypoints)

    expect_true(is.numeric(result))
    expect_true(length(result) > 0)
    expect_true(all(result >= 1 & result <= nrow(wp_data)))
})

# Additional tests for edge cases
test_that(".terminalStateFromMarkovChain handles edge cases", {
    # Test with minimum waypoints
    waypoints_min <- sample(seq_len(50), 2)
    Transmat_min <- matrix(runif(4), nrow = 2, ncol = 2)
    Transmat_min <- Transmat_min / rowSums(Transmat_min)
    result_min <- .terminalStateFromMarkovChain(Transmat_min, wp_data, pseudotime, waypoints_min)
    expect_true(is.numeric(result_min))
    expect_true(length(result_min) > 0)

    # Test with maximum waypoints
    waypoints_max <- sample(seq_len(50), 50)
    Transmat_max <- matrix(runif(2500), nrow = 50, ncol = 50)
    Transmat_max <- Transmat_max / rowSums(Transmat_max)
    result_max <- .terminalStateFromMarkovChain(Transmat_max, wp_data, pseudotime, waypoints_max)
    expect_true(is.numeric(result_max))
    expect_true(length(result_max) > 0)
})
