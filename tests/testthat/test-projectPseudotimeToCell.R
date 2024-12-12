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
# Construct Pseudobulked VDJ Feature Space
pb.milo <- vdjPseudobulk(milo_object, col_to_take = "anno_lvl_2_final_clean")
pb.milo <- scater::runPCA(pb.milo, assay.type = "Feature_space")
# Define root and branch tips
pca <- t(as.matrix(SingleCellExperiment::reducedDim(pb.milo, type = "PCA")))
branch.tips <- c(189, 198) # which.min(pca[, 2]) and which.max(pca[, 2])
names(branch.tips) <- c("CD8+T", "CD4+T")
root <- 177 # which.max(pca[, 1])
# Construct Diffusion Map
dm <- destiny::DiffusionMap(t(pca), n_pcs = 50, n_eigs = 10)
dif.pse <- destiny::DPT(dm, tips = c(root, branch.tips), w_width = 0.1)
pb.milo <- markovProbability(
    milo = pb.milo,
    diffusionmap = dm,
    terminal_state = branch.tips,
    diffusiontime = dif.pse[[paste0("DPT", root)]],
    root_cell = root,
    pseudotime_key = "pseudotime"
)
# Test for projectPseudotimeToCell
test_that("projectPseudotimeToCell works correctly", {
    result <- projectPseudotimeToCell(milo_object, pb.milo, branch.tips)

    expect_true(is(result, "SingleCellExperiment"))
    expect_true("pseudotime" %in% colnames(colData(result)))
    expect_true(all(c("CD8+T", "CD4+T") %in% colnames(colData(result))))
})

# Additional tests for edge cases
test_that("projectPseudotimeToCell cannot work without pseudotime column", {
    # Test with wrong column
    expect_error(projectPseudotimeToCell(milo_object, pb.milo, branch.tips, pseudotime_key = "wrong_pseudotime"), "subscript contains invalid names")
    # Test with missing pseudotime
    colData(pb.milo)$pseudotime <- NULL
    expect_error(projectPseudotimeToCell(milo_object, pb.milo, branch.tips), "subscript contains invalid names")
})
