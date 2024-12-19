library(scRepertoire)
library(miloR)
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
sce <- setupVdjPseudobulk(sce,
    mode_option = "abT",
    already.productive = TRUE,
    subsetby = "anno_lvl_2_final_clean",
    groups = c("CD8+T", "CD4+T", "ABT(ENTRY)", "DP(P)_T", "DP(Q)_T")
)

# Build Milo Object
milo_object <- Milo(sce)
milo_object <- buildGraph(milo_object, k = 30, d = 20, reduced.dim = "X_scvi")
milo_object <- makeNhoods(milo_object, reduced_dims = "X_scvi", d = 20, prop = 0.3)

# Test for miloUmap
test_that("miloUmap works correctly", {
    result <- miloUmap(milo_object, n_neighbors = 30)

    expect_true(is(result, "Milo"))
    expect_true("UMAP_knngraph" %in% reducedDimNames(result))
    expect_equal(ncol(reducedDim(result, "UMAP_knngraph")), 2)
})
