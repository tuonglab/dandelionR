library(destiny)

# Mock data for testing
set.seed(100)
data <- matrix(rnorm(100), ncol = 10)
dm <- DiffusionMap(data)

# Test for .determineMultiscaleSpace
test_that(".determineMultiscaleSpace works correctly", {
    result <- .determineMultiscaleSpace(dm)

    expect_true(is.data.frame(result))
    expect_true(ncol(result) > 0)
    expect_true(nrow(result) == nrow(data))
})

# Test with specified n_eigs
test_that(".determineMultiscaleSpace works with specified n_eigs", {
    result <- .determineMultiscaleSpace(dm, n_eigs = 5)

    expect_true(is.data.frame(result))
    expect_true(ncol(result) == 5)
    expect_true(nrow(result) == nrow(data))
})

# Additional tests for edge cases
test_that(".determineMultiscaleSpace handles edge cases", {
    # Test with minimum n_eigs
    result_min <- .determineMultiscaleSpace(dm, n_eigs = 1)
    expect_true(is.data.frame(result_min))
    expect_true(ncol(result_min) == 1)
    expect_true(nrow(result_min) == nrow(data))

    # Test with maximum n_eigs
    result_max <- .determineMultiscaleSpace(dm, n_eigs = ncol(data))
    expect_true(is.data.frame(result_max))
    expect_true(ncol(result_max) == ncol(data))
    expect_true(nrow(result_max) == nrow(data))
})
