library(SingleCellExperiment)

# Mock data for testing
set.seed(100)
data(sce_vdj)

# Test for setupVdjPseudobulk
test_that("setupVdjPseudobulk works correctly with mode_option", {
    result <- setupVdjPseudobulk(
        sce = sce_vdj,
        mode_option = "abT",
        already.productive = FALSE,
        productive_cols = NULL,
        productive_vj = TRUE,
        productive_vdj = TRUE,
        allowed_chain_status = c("Single pair", "Extra pair")
    )

    expect_true(is(result, "SingleCellExperiment"))
    expect_true(ncol(result) <= ncol(sce_vdj))
})

test_that("setupVdjPseudobulk works correctly with productive_cols", {
    result <- setupVdjPseudobulk(
        sce = sce_vdj,
        mode_option = NULL,
        already.productive = FALSE,
        productive_cols = c("productive_abT_VDJ", "productive_abT_VJ"),
        productive_vj = TRUE,
        productive_vdj = TRUE,
        allowed_chain_status = c("Single pair", "Extra pair")
    )

    expect_true(is(result, "SingleCellExperiment"))
    expect_true(ncol(result) <= ncol(sce_vdj))
})

test_that("setupVdjPseudobulk handles allowed_chain_status correctly", {
    result <- setupVdjPseudobulk(
        sce = sce_vdj,
        mode_option = "abT",
        already.productive = TRUE,
        allowed_chain_status = c("Single pair")
    )

    expect_true(is(result, "SingleCellExperiment"))
    expect_true(all(colData(result)$chain_status == "Single pair"))
})

test_that(sprintf(
    "setupVdjPseudobulk throws error when no cells match %s",
    "allowed_chain_status"
), {
    expect_error(setupVdjPseudobulk(
        sce = sce_vdj,
        mode_option = "abT",
        already.productive = TRUE,
        allowed_chain_status = c("Non-existent status")
    ), "Unsuitable allowed_chain_status")
})
