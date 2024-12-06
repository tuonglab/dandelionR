#' projectPseudotimeToCell
#'
#' Function to project pseudotime & branch probabilities from pb.milo (pseudobulk) to milo (cell).
#' @param milo SingleCellExperiment or milo object
#' @param pb_milo pseudobulk data
#' @param term_states vector of terminal states with branch_probabilities to be transferred
#' @param suffix suffix to be added after the added column names, default ''
#' @examples
#' data(sce_vdj)
#' sce_vdj <- setupVdjPseudobulk(sce_vdj,
#'     already.productive = FALSE
#' )
#' # Build Milo Object
#' set.seed(100)
#' traj_milo <- miloR::Milo(sce_vdj)
#' milo_object <- miloR::buildGraph(traj_milo, k = 50, d = 20, reduced.dim = "X_scvi")
#' milo_object <- miloR::makeNhoods(milo_object, reduced_dims = "X_scvi", d = 20)
#'
#' # Construct Pseudobulked VDJ Feature Space
#' pb.milo <- vdjPseudobulk(milo_object, col_to_take = "anno_lvl_2_final_clean")
#' pb.milo <- scater::runPCA(pb.milo, assay.type = "Feature_space")
#'
#' # Define root and branch tips
#' pca <- t(as.matrix(SingleCellExperiment::reducedDim(pb.milo, type = "PCA")))
#' branch.tips <- c(232, 298)
#' names(branch.tips) <- c("CD8+T", "CD4+T")
#' root <- 476
#'
#' # Construct Diffusion Map
#' dm <- destiny::DiffusionMap(t(pca), n_pcs = 50, n_eigs = 10)
#' dif.pse <- destiny::DPT(dm, tips = c(root, branch.tips), w_width = 0.1)
#'
#' # Markov Chain Construction
#' pb.milo <- markovProbability(
#'     milo = pb.milo,
#'     diffusionmap = dm,
#'     terminal_state = branch.tips,
#'     diffusiontime = dif.pse[[paste0("DPT", root)]],
#'     root_cell = root,
#'     pseudotime_key = "pseudotime"
#' )
#'
#' # Project Pseudobulk Data
#' projected_milo <- projectPseudotimeToCell(milo_object, pb.milo, branch.tips)
#'
#' @return subset of milo or SingleCellExperiment object where cell that do not belong to any neighbourhood are removed and projected pseudotime information stored colData
#' @import miloR
#' @import SingleCellExperiment
#' @export
projectPseudotimeToCell <- function(milo, pb_milo, term_states, suffix = "") {
    nhood <- nhoods(pb_milo) # peudobulk x cells
    # leave out cells that do not blongs to any neighbourhood
    nhoodsum <- apply(nhoods(pb_milo), 2, sum)
    cdata <- milo[, nhoodsum > 0]
    message(sprintf(
        "%d number of cells removed due to not belonging to any neighbourhood",
        sum(nhoodsum == 0)
    ))
    # for each cell pesudotime_mean is the average of the pseudobulks the cell
    # is in, weighted by 1/ neighbourhood size
    nhoods_cdata <- nhood[, nhoodsum > 0]
    nhoods_cdata_norm <- nhoods_cdata / apply(nhoods_cdata, 1, sum)
    col_list <- c(paste0("pseudotime", suffix), paste0(names(term_states), suffix))
    new_col <- vapply(colData(pb_milo)[, col_list]@listData, function(x, y) {
        list(as.vector(x %*% y / apply(y, 2, sum)))
    }, y = nhoods_cdata_norm, FUN.VALUE = list(double()))
    colData(cdata) <- cbind(colData(cdata), new_col)
    cdata
}
