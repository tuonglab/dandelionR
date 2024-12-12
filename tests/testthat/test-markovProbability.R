library(miloR)
library(destiny)
library(SingleCellExperiment)

# Mock data for testing
set.seed(100)
data <- matrix(rnorm(100), ncol = 10)
sce <- SingleCellExperiment(assays = list(counts = data))
colData(sce)$pseudotime <- runif(ncol(sce))

# Build Milo Object
milo_object <- Milo(sce)
milo_object <- buildGraph(milo_object, k = 5, d = 2, reduced.dim = "counts")
milo_object <- makeNhoods(milo_object, reduced_dims = "counts", d = 2)

# Construct Pseudobulked VDJ Feature Space
pb.milo <- vdjPseudobulk(milo_object, col_to_take = "pseudotime")
pb.milo <- scater::runPCA(pb.milo, assay.type = "counts")

# Define root and branch tips
pca <- t(as.matrix(SingleCellExperiment::reducedDim(pb.milo, type = "PCA")))
branch.tips <- c(which.min(pca[2, ]), which.max(pca[2, ]))
names(branch.tips) <- c("CD8+T", "CD4+T")
root <- which.max(pca[1, ])

# Construct Diffusion Map
dm <- DiffusionMap(t(pca), n_pcs = 5, n_eigs = 5)
dif.pse <- destiny::DPT(dm, tips = c(root, branch.tips), w_width = 0.1)

# Test for markovProbability
test_that("markovProbability works correctly", {
    result <- markovProbability(
        milo = pb.milo,
        diffusionmap = dm,
        diffusiontime = dif.pse[[paste0("DPT", root)]],
        terminal_state = branch.tips,
        root_cell = root,
        pseudotime_key = "pseudotime"
    )

    expect_true(is(result, "SingleCellExperiment"))
    expect_true("pseudotime" %in% colnames(colData(result)))
    expect_true(all(c("CD8+T", "CD4+T") %in% colnames(colData(result))))
})

# Additional tests for edge cases
test_that("markovProbability handles edge cases", {
    # Test with missing pseudotime
    colData(pb.milo)$pseudotime <- NULL
    result_no_pseudotime <- markovProbability(
        milo = pb.milo,
        diffusionmap = dm,
        diffusiontime = dif.pse[[paste0("DPT", root)]],
        terminal_state = branch.tips,
        root_cell = root,
        pseudotime_key = "pseudotime"
    )
    expect_true(is(result_no_pseudotime, "SingleCellExperiment"))
    expect_true("pseudotime" %in% colnames(colData(result_no_pseudotime)))

    # Test with minimum waypoints
    result_min_waypoints <- markovProbability(
        milo = pb.milo,
        diffusionmap = dm,
        diffusiontime = dif.pse[[paste0("DPT", root)]],
        terminal_state = branch.tips,
        root_cell = root,
        pseudotime_key = "pseudotime",
        num_waypoints = 2
    )
    expect_true(is(result_min_waypoints, "SingleCellExperiment"))

    # Test with maximum waypoints
    result_max_waypoints <- markovProbability(
        milo = pb.milo,
        diffusionmap = dm,
        diffusiontime = dif.pse[[paste0("DPT", root)]],
        terminal_state = branch.tips,
        root_cell = root,
        pseudotime_key = "pseudotime",
        num_waypoints = 50
    )
    expect_true(is(result_max_waypoints, "SingleCellExperiment"))
})
