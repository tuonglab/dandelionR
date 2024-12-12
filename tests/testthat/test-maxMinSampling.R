# Mock data for testing
set.seed(123)
data <- matrix(rnorm(100), ncol = 5)

# Test for .maxMinSampling
test_that(".maxMinSampling works correctly", {
    num_waypoints <- 10
    result <- .maxMinSampling(data, num_waypoints)

    expect_true(is.numeric(result))
    expect_true(length(result) <= num_waypoints)
    expect_true(all(result >= 1 & result <= nrow(data)))
})

# Additional tests for edge cases
test_that(".maxMinSampling handles edge cases", {
    # Test with minimum waypoints
    result_min <- .maxMinSampling(data, 1)
    expect_true(is.numeric(result_min))
    expect_true(length(result_min) == 1)

    # Test with maximum waypoints
    result_max <- .maxMinSampling(data, nrow(data))
    expect_true(is.numeric(result_max))
    expect_true(length(result_max) <= nrow(data))

    # Test with single column data
    single_col_data <- matrix(rnorm(20), ncol = 1)
    result_single_col <- .maxMinSampling(single_col_data, 5)
    expect_true(is.numeric(result_single_col))
    expect_true(length(result_single_col) <= 5)

    # Test with single row data
    single_row_data <- matrix(rnorm(5), nrow = 1)
    result_single_row <- .maxMinSampling(single_row_data, 1)
    expect_true(is.numeric(result_single_row))
    expect_true(length(result_single_row) == 1)
})
