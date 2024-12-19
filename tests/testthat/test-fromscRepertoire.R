library(SingleCellExperiment)

# Mock data for testing
library(dandelionR)
library(scRepertoire)

data(demo_sce)
data(demo_airr)

contig.list <- loadContigs(input = demo_airr, format = "AIRR")

# Format to `scRepertoire`'s requirements and some light filtering
combined.TCR <- combineTCR(contig.list,
    removeNA = TRUE,
    removeMulti = FALSE,
    filterMulti = TRUE
)
sce <- combineExpression(combined.TCR, demo_sce)

# Test for setupVdjPseudobulk
test_that("setupVdjPseudobulk works correctly", {
    result <- setupVdjPseudobulk(sce,
        mode_option = "abT",
        already.productive = TRUE,
        subsetby = "anno_lvl_2_final_clean",
        groups = c("CD8+T", "CD4+T", "ABT(ENTRY)", "DP(P)_T", "DP(Q)_T")
    )

    expect_true(is(result, "SingleCellExperiment"))
    expect_true(ncol(result) <= ncol(sce))
    expect_true(all(colData(result)$anno_lvl_2_final_clean %in% c("CD8+T", "CD4+T", "ABT(ENTRY)", "DP(P)_T", "DP(Q)_T")))
})
