library(SingleCellExperiment)

# Mock data for testing
set.seed(100)
data <- matrix(rnorm(100), ncol = 10)
col_data <- data.frame(
    CTgene = c(
        "TRA1.TRA2_TRB1.TRB2.TRB3", "TRA1.TRA2_TRB1.TRB2.TRB3", "NA_NA",
        "TRA1.TRA2_TRB1.TRB2.TRB3", "TRA1.TRA2_TRB1.TRB2.TRB3", "NA_NA",
        "TRA1.TRA2_TRB1.TRB2.TRB3", "TRA1.TRA2_TRB1.TRB2.TRB3", "NA_NA",
        "TRA1.TRA2_TRB1.TRB2.TRB3"
    )
)
sce <- SingleCellExperiment(assays = list(counts = data), colData = col_data)

# Test for splitCTgene
test_that("splitCTgene works correctly", {
    result <- splitCTgene(sce)

    expect_true(is.list(result))
    expect_equal(length(result), ncol(sce))
    expect_true(all(sapply(result, length) == 5))
})

# Test for formatVdj
test_that("formatVdj works correctly", {
    gene_list <- list(c("TRA1", "TRA2"), c("TRB1", "TRB2", "TRB3"))
    result <- formatVdj(gene_list)

    expect_true(is.character(result))
    expect_equal(length(result), 5)
    expect_equal(result, c("TRA1", "TRA2", "TRB1", "TRB2", "TRB3"))
})

# Test for chainAssign
test_that("chainAssign works correctly", {
    vec <- c("TRA1", "TRA2", "TRB1", "TRB2", "TRB3")
    result_vj <- chainAssign(vec[1:2], 2)
    result_vdj <- chainAssign(vec[3:5], 3)

    expect_true(is.character(result_vj))
    expect_equal(length(result_vj), 2)
    expect_equal(result_vj, c("TRA1", "TRA2"))

    expect_true(is.character(result_vdj))
    expect_equal(length(result_vdj), 3)
    expect_equal(result_vdj, c("TRB1", "TRB2", "TRB3"))
})

# Additional tests for edge cases
test_that("splitCTgene handles edge cases", {
    colData(sce)$CTgene <- c("NA_NA", "NA_NA", "NA_NA", "NA_NA", "NA_NA", "NA_NA", "NA_NA", "NA_NA", "NA_NA", "NA_NA")
    result_na <- splitCTgene(sce)
    expect_true(is.list(result_na))
    expect_equal(length(result_na), ncol(sce))
    expect_true(all(sapply(result_na, function(x) all(x == "None"))))
})

test_that("chainAssign handles edge cases", {
    vec_na <- c("NA", "NA")
    result_na <- chainAssign(vec_na, 2)
    expect_true(is.character(result_na))
    expect_equal(length(result_na), 2)
    expect_equal(result_na, c("None", "None"))

    expect_error(chainAssign(vec_na, 3), "argument \"num\" is missing, with no default")
})
