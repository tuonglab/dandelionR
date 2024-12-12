# Mock data for testing
set.seed(100)
data <- matrix(abs(rnorm(100)), ncol = 10)

# Test for .minMaxScale
test_that(".minMaxScale works correctly", {
    result <- .minMaxScale(data)

    expect_true(is.matrix(result))
    expect_equal(dim(result), dim(data))
    expect_true(all(result >= 0 & result <= 1))

    # Check if the minimum value in each column is 0 and the maximum value is 1
    expect_equal(apply(result, 2, min), rep(0, ncol(data)))
    expect_equal(apply(result, 2, max), rep(1, ncol(data)))
})

# Additional tests for edge cases
test_that(".minMaxScale handles edge cases", {
    # Test with single column data
    single_col_data <- matrix(abs(rnorm(10)), ncol = 1)
    result_single_col <- .minMaxScale(single_col_data)
    expect_true(is.matrix(result_single_col))
    expect_equal(dim(result_single_col), rev(dim(single_col_data)))
    expect_true(all(result_single_col >= 0 & result_single_col <= 1))

    # Test with single row data
    single_row_data <- matrix(abs(rnorm(10)), nrow = 1)
    result_single_row <- .minMaxScale(single_row_data)
    expect_true(all(is.na(result_single_row)))

    # Test with constant data
    constant_data <- matrix(rep(5, 100), ncol = 10)
    result_constant <- .minMaxScale(constant_data)
    expect_true(is.matrix(result_constant))
    expect_equal(dim(result_constant), dim(constant_data))
    expect_true(all(is.na(result_constant)))
})
