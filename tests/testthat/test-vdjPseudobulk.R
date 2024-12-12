library(miloR)
library(SingleCellExperiment)

# Mock data for testing
set.seed(100)
data(sce_vdj)
sce_vdj <- setupVdjPseudobulk(sce_vdj,
    already.productive = FALSE,
    allowed_chain_status = c("Single pair", "Extra pair")
)
# Build Milo Object
milo_object <- miloR::Milo(sce_vdj)
milo_object <- miloR::buildGraph(milo_object, k = 50, d = 20, reduced.dim = "X_scvi")
milo_object <- miloR::makeNhoods(milo_object, reduced_dims = "X_scvi", d = 20)

# Test for vdjPseudobulk
test_that("vdjPseudobulk works correctly with Milo object", {
    result <- vdjPseudobulk(milo_object, col_to_take = "anno_lvl_2_final_clean")

    expect_true(is(result, "SingleCellExperiment"))
    expect_true("Feature_space" %in% assayNames(result))
    expect_true(ncol(result) > 0)
    expect_true(nrow(result) > 0)
})

test_that("vdjPseudobulk doesn't work with just any sce object", {
    expect_error(vdjPseudobulk(sce_vdj, col_to_bulk = "anno_lvl_2_final_clean"), "object 'pbs.col' not found")
})
