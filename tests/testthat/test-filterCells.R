library(rlang)

# Mock data for testing
set.seed(100)
data <- matrix(rnorm(100), ncol = 10)
col_data <- data.frame(
    VDJ = c(
        "None", "TRBV7-9", "TRBV7-9", "None", "TRBV7-9",
        "No_contig", "TRBV7-9", "TRBV7-9", "None", "TRBV7-9"
    )
)
sce <- SingleCellExperiment(assays = list(counts = data), colData = col_data)

# Test for .filterCells
test_that(".filterCells works correctly with remove_missing = TRUE", {
    result <- .filterCells(sce, "VDJ", remove_missing = TRUE)

    expect_true(is(result, "SingleCellExperiment"))
    expect_equal(ncol(result), 6)
    expect_false(any(colData(result)$VDJ %in% c("None", "No_contig")))
})

test_that(".filterCells works correctly with remove_missing = FALSE", {
    result <- .filterCells(sce, "VDJ", remove_missing = FALSE)

    expect_true(is(result, "SingleCellExperiment"))
    expect_equal(ncol(result), 10)
    expect_true(all(colData(result)$VDJ[grep(
        "None|No_contig",
        colData(result)$VDJ
    )] == "VDJ_missing"))
})

# Additional tests for edge cases
test_that(".filterCells handles edge cases", {
    # Test with empty SingleCellExperiment
    empty_sce <- SingleCellExperiment(assays = list(counts = matrix(
        nrow = 0,
        ncol = 0
    )))
    expect_error(
        .filterCells(empty_sce, "VDJ"),
        sprintf(paste(
            "None column remains, please check whether the",
            "filtering option is correct."
        ))
    )

    # Test with no matching filter pattern
    no_match_sce <- SingleCellExperiment(
        assays = list(counts = data),
        colData = data.frame(VDJ = rep("TRBV7-9", 10))
    )
    result_no_match <- .filterCells(no_match_sce, "VDJ", remove_missing = TRUE)
    expect_equal(ncol(result_no_match), 10)
})
