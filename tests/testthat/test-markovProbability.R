# Mock data for testing
data(sce_vdj)
# downsample to first 2000 cells to speed up the testing.
sce_vdj <- sce_vdj[, 1:2000]
sce_vdj <- setupVdjPseudobulk(sce_vdj,
    already.productive = FALSE,
    allowed_chain_status = c("Single pair", "Extra pair")
)

# Build Milo Object
set.seed(100)
milo_object <- miloR::Milo(sce_vdj)
milo_object <- miloR::buildGraph(milo_object,
    k = 20, d = 20,
    reduced.dim = "X_scvi"
)
milo_object <- miloR::makeNhoods(milo_object, reduced_dims = "X_scvi", d = 20)

# Construct Pseudobulked VDJ Feature Space
pb.milo <- vdjPseudobulk(milo_object, col_to_take = "anno_lvl_2_final_clean")
pb.milo <- scater::runPCA(pb.milo, assay.type = "Feature_space")

# Define root and branch tips
pca <- t(as.matrix(SingleCellExperiment::reducedDim(pb.milo, type = "PCA")))
branch.tips <- c(which.min(pca[, 2]), which.max(pca[, 2]))
names(branch.tips) <- c("CD8+T", "CD4+T")
root <- which.min(pca[, 1])

# Construct Diffusion Map
dm <- destiny::DiffusionMap(t(pca), n_pcs = 10, n_eigs = 5)
dif.pse <- destiny::DPT(dm, tips = c(root, branch.tips), w_width = 0.1)
# Since the markovProbability will also run in test-projectPseudotimeToCell.R
# as setup, we only consider two special circumstances here.
test_that("markovProbability works correctly with terminal states as NULL", {
    pb.milo <- markovProbability(
        milo = pb.milo,
        diffusionmap = dm,
        terminal_state = NULL,
        diffusiontime = dif.pse[[paste0("DPT", root)]],
        root_cell = root,
        pseudotime_key = "pseudotime"
    )
    expect_true(is(pb.milo, "Milo"))
    expect_true("pseudotime" %in% colnames(colData(pb.milo)))
    expect_true(any(grep("terminal_state", colnames(colData(pb.milo)))))
    expect_false(is.null(metadata(pb.milo)))
})

test_that("markovProbability works correctly with only one terminal state", {
    pb.milo <- markovProbability(
        milo = pb.milo,
        diffusionmap = dm,
        terminal_state = as.numeric(branch.tips[1]),
        diffusiontime = dif.pse[[paste0("DPT", root)]],
        root_cell = root,
        pseudotime_key = "pseudotime"
    )
    expect_true(is(pb.milo, "Milo"))
    expect_true("pseudotime" %in% colnames(colData(pb.milo)))
    expect_true(any(grep("terminal_state", colnames(colData(pb.milo)))))
    expect_false(is.null(metadata(pb.milo)))
})
