test_that("test the full workflow", {
    data(sce_vdj)
    sce_vdj <- setupVdjPseudobulk(sce_vdj,
        already.productive = FALSE
    )
    # Build Milo Object
    set.seed(100)
    traj_milo <- miloR::Milo(sce_vdj)
    milo_object <- miloR::buildGraph(traj_milo, k = 50, d = 20, reduced.dim = "X_scvi")
    milo_object <- miloR::makeNhoods(milo_object, reduced_dims = "X_scvi", d = 20)

    # Construct Pseudobulked VDJ Feature Space
    pb.milo <- vdjPseudobulk(milo_object, col_to_take = "anno_lvl_2_final_clean")
    pb.milo <- scater::runPCA(pb.milo, assay.type = "Feature_space")

    # Define root and branch tips
    pca <- t(as.matrix(SingleCellExperiment::reducedDim(pb.milo, type = "PCA")))
    branch.tips <- c(232, 298)
    names(branch.tips) <- c("CD8+T", "CD4+T")
    root <- 476

    # Construct Diffusion Map
    dm <- destiny::DiffusionMap(t(pca), n_pcs = 50, n_eigs = 10)
    dif.pse <- destiny::DPT(dm, tips = c(root, branch.tips), w_width = 0.1)

    # Markov Chain Construction
    pb.milo <- markovProbability(
        milo = pb.milo,
        diffusionmap = dm,
        terminal_state = branch.tips,
        diffusiontime = dif.pse[[paste0("DPT", root)]],
        root_cell = root,
        pseudotime_key = "pseudotime"
    )

    # Project Pseudobulk Data
    projected_milo <- projectPseudotimeToCell(milo_object, pb.milo, branch.tips)
})
