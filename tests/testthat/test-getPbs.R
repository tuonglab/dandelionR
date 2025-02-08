library(SingleCellExperiment)
library(Matrix)
library(rlang)

# Mock data for testing
set.seed(100)
data <- matrix(rnorm(100), ncol = 10)
col_data <- data.frame(
    sample_type = rep(c("A", "B"), each = 5),
    condition = rep(c("control", "treated"), times = 5)
)
sce <- SingleCellExperiment(assays = list(counts = data), colData = col_data)

# Test for .getPbs
test_that(".getPbs works correctly with pbs provided", {
    pbs <- Matrix::Matrix(matrix(1, nrow = 10, ncol = 2), sparse = TRUE)
    result <- .getPbs(pbs, NULL, sce)

    expect_true(is(result, "Matrix"))
    expect_equal(dim(result), c(10, 2))
})

test_that(".getPbs works correctly with col_to_bulk provided", {
    result <- .getPbs(NULL, c("sample_type", "condition"), sce)

    expect_true(is(result, "Matrix"))
    expect_equal(dim(result), c(10, 4))
    expect_equal(colnames(result), c(
        "A,control", "A,treated",
        "B,control", "B,treated"
    ))
})

test_that(".getPbs throws error when both pbs and col_to_bulk are NULL", {
    expect_error(
        .getPbs(NULL, NULL, sce),
        "You must specify 'pbs' or 'col_to_bulk'."
    )
})

test_that(".getPbs throws error when both pbs and col_to_bulk are provided", {
    pbs <- Matrix::Matrix(matrix(1, nrow = 10, ncol = 2), sparse = TRUE)
    expect_error(
        .getPbs(pbs, c("sample_type", "condition"), sce),
        "You must specify 'pbs' or 'col_to_bulk', not both."
    )
})
