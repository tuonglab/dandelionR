#' markovProbability
#'
#' Preprocessing data and Construct markov chain and calculate probabilities
#' @param milo milo or SingelCellExperiment object, with pseudotime stored in colData, used to store the result and extract pseudotime. Pseudotime stored in milo has higher priority than the value provided through the diffusiontime parameter
#' @param diffusionmap DiffusionMap object corresponds to milo
#' @param diffusiontime if milo do not restore pseudotime, use this parameter to transfer it to function
#' @param terminal_state the index of the terminal state
#' @param root_cell the index of the root state
#' @param pseudotime_key the column name in the colData that holds the inferred pseudotime
#' @param scale_components logical, If True, the components will be scale before constructing markov chain
#' @param num_waypoints integer, 500L by default. Number of waypoints to sample to construct markov chain.
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
#'     diffusiontime = dif.pse[[paste0("DPT", root)]],
#'     terminal_state = branch.tips,
#'     root_cell = root,
#'     pseudotime_key = "pseudotime"
#' )
#'
#' @return milo or SinglCellExperiment object with pseudotime, probabilities in its colData
#' @include determMultiscaleSpace.R
#' @include minMaxScale.R
#' @include maxMinSampling.R
#' @include differentiationProbabilities.R
#' @include projectProbability.R
#' @import SingleCellExperiment
#' @export
markovProbability <- function(
    milo, diffusionmap, diffusiontime = NULL, terminal_state, root_cell,
    pseudotime_key = "pseudotime",
    scale_components = TRUE, num_waypoints = 500) {
    if (is.null(milo[[pseudotime_key]])) {
        requireNamespace("rlang")
        if (is.null(diffusiontime)) {
            rlang::abort(paste("Missing pseudotime data. This data can be either stored in", deparse(substitute(milo)), "or provided by parameter diffusiontime"))
        } else {
            milo[[pseudotime_key]] <- diffusiontime
        }
    } else {
        diffusiontime <- milo[[pseudotime_key]]
    }

    # scale data
    multiscale <- .determineMultiscaleSpace(diffusionmap)
    if (scale_components) {
        multiscale <- .minMaxScale(multiscale)
    }
    # sample waypoints to construct markov chain
    waypoints <- .maxMinSampling(multiscale, num_waypoints = 500)
    waypoints <- unique(c(root_cell, waypoints, terminal_state))
    # calculate probabilities
    probabilities <- differentiationProbabilities(multiscale[waypoints, ],
        terminal_states = terminal_state,
        pseudotime = diffusiontime, waypoints = waypoints
    )
    # project probabilities from waypoints to each pseudobulk
    probabilities_proj <- projectProbability(diffusionmap, waypoints, probabilities)
    # store the result into milo
    requireNamespace("S4Vectors")
    new_coldata <- S4Vectors::DataFrame(probabilities_proj[
        ,
        1
    ], probabilities_proj[, 2])
    colnames(new_coldata) <- c(names(terminal_state))
    # prevent same name in colData
    idx <- names(colData(milo)) %in% colnames(new_coldata)
    if (any(idx)) {
        warning(paste(
            "Name", paste(names(colData(milo))[idx], collapse = ", "),
            "already exists in", as.character(substitute(milo))
        ))
        repeat {
            answer <- readline(prompt = "Do you want to overwrite the column? (y/n): ")
            if (answer == "n") {
                while (any(names(colData(milo)) %in% colnames(new_coldata))) {
                    colnames(new_coldata) <- paste0(colnames(new_coldata), "_new")
                }
                msg <- paste(colnames(new_coldata), collapse = ", ")
                message(sprintf("The data will stored in %s", msg))
                break
            } else if (answer == "y") {
                msg <- paste(names(colData(milo))[idx], collapse = ", ")
                message(sprintf("Overwriting %s ...", msg))
                colData(milo) <- colData(milo)[!idx]
                break
            } else {
                message("Invalid response. Please enter 'y' or 'n'.")
            }
        }
    }
    colData(milo) <- cbind(colData(milo), new_coldata)
    return(milo)
}
