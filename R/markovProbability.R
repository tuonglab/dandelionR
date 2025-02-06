#' Markov Chain Construction and Probability Calculation
#'
#' This function preprocesses data, constructs a Markov chain, and calculates transition probabilities based on pseudotime information.
#' @param milo A `Milo` or `SingleCellExperiment` object. This object should have pseudotime stored in `colData`, which will be used to calculate probabilities. If pseudotime is available in `milo`, it takes precedence over the value provided through the `diffusiontime` parameter.
#' @param terminal_state Integer. The index of the terminal state in the Markov chain.
#' @param root_cell Integer. The index of the root state in the Markov chain.
#' @param knn Integer. The number of nearest neighbors for graph construction. Default is `30L`.
#' @param diffusionmap A `DiffusionMap` object corresponding to the `milo` object. Used for Markov chain construction.
#' @param diffusiontime Numeric vector. If pseudotime is not stored in `milo`, this parameter can be used to provide pseudotime values to the function.
#' @param pseudotime_key Character. The name of the column in `colData` that contains the inferred pseudotime.
#' @param scale_components Logical. If `TRUE`, the components will be scaled before constructing the Markov chain. Default is `FALSE`.
#' @param num_waypoints Integer. The number of waypoints to sample when constructing the Markov chain. Default is `500L`.
#' @param verbose Boolean, whether to print messages/warnings.
#' @examples
#' data(sce_vdj)
#' sce_vdj <- setupVdjPseudobulk(sce_vdj,
#'     already.productive = FALSE,
#'     allowed_chain_status = c("Single pair", "Extra pair")
#' )
#' # Build Milo Object
#' set.seed(100)
#' milo_object <- miloR::Milo(sce_vdj)
#' milo_object <- miloR::buildGraph(milo_object, k = 50, d = 20, reduced.dim = "X_scvi")
#' milo_object <- miloR::makeNhoods(milo_object, reduced_dims = "X_scvi", d = 20)
#'
#' # Construct Pseudobulked VDJ Feature Space
#' pb.milo <- vdjPseudobulk(milo_object, col_to_take = "anno_lvl_2_final_clean")
#' pb.milo <- scater::runPCA(pb.milo, assay.type = "Feature_space")
#'
#' # Define root and branch tips
#' pca <- t(as.matrix(SingleCellExperiment::reducedDim(pb.milo, type = "PCA")))
#' branch.tips <- c(189, 198) # which.min(pca[, 2]) and which.max(pca[, 2])
#' names(branch.tips) <- c("CD8+T", "CD4+T")
#' root <- 177 # which.max(pca[, 1])
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
#' @importFrom rlang abort
#' @importFrom S4Vectors DataFrame metadata metadata<-
#' @importFrom SummarizedExperiment colData<-
#' @export
markovProbability <- function(
    milo, diffusionmap, terminal_state = NULL, root_cell, knn = 30L,
    diffusiontime = NULL, pseudotime_key = "pseudotime", scale_components = TRUE,
    num_waypoints = 500, verbose = TRUE) {
    if (is.null(milo[[pseudotime_key]])) {
        if (is.null(diffusiontime)) { # nocov start
            abort(paste(
                "Missing pseudotime data. This data can be either stored in",
                deparse(substitute(milo)), "or provided by parameter diffusiontime"
            )) # nocov end
        } else {
            milo[[pseudotime_key]] <- diffusiontime
        }
    } else {
        diffusiontime <- milo[[pseudotime_key]] # nocov
    }

    # scale data
    multiscale <- .determineMultiscaleSpace(diffusionmap)
    if (scale_components) {
        multiscale <- .minMaxScale(multiscale)
    }
    # sample waypoints to construct markov chain
    waypoints <- .maxMinSampling(multiscale, num_waypoints = 500, verbose = verbose)
    waypoints <- unique(c(root_cell, waypoints, terminal_state))
    # calculate probabilities
    probabilities_terminal <- differentiationProbabilities(multiscale[waypoints, ],
        terminal_states = terminal_state,
        knn = knn, pseudotime = diffusiontime, waypoints = waypoints, verbose = verbose
    )
    probabilities <- probabilities_terminal[[1]]
    terminal_state <- probabilities_terminal[[2]]
    # project probabilities from waypoints to each pseudobulk
    probabilities_proj <- projectProbability(diffusionmap, waypoints, probabilities, verbose)
    # store the result into milo
    milo <- .addColData(probabilities_proj, terminal_state, milo)
    return(milo)
}


.addColData <- function(probabilities_proj, terminal_state, milo){
  new_coldata <- DataFrame(as.matrix(probabilities_proj))
  if (is.null(names(terminal_state))) {
    names(terminal_state) <- paste0("terminal_state", seq(length(terminal_state)))
  }
  colnames(new_coldata) <- c(names(terminal_state))
  # prevent same name in colData
  idx <- names(colData(milo)) %in% colnames(new_coldata)
  if (any(idx)) { # nocov start
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
  } # nocov end
  colData(milo) <- cbind(colData(milo), new_coldata)
  metadata(milo) <- c(metadata(milo), list(branch.tips = terminal_state))
  return(milo)
}
