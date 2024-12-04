#' markov_probability
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
#' @return milo or SinglCellExperiment object with pseudotime, probabilities in its colData
#' @include determ.multiscale.space.R
#' @include minmax.scale.R
#' @include max.min.sampling.R
#' @include differentiation_probabilities.R
#' @include project_probability.R
#' @import SingleCellExperiment
#' @examples
#' # load data
#' data(pb.milo)
#' library(SingleCellExperiment)
#' branch.tips <- c(540,54)
#' names(branch.tips) <- c("CD8+T", "CD4+T")
#' root <- 222
#' DPTroot <- paste0("DPT", root)
#' # extract pca
#' pca <- t(as.matrix(reducedDim(pb.milo, type = "PCA")))
#' # Run diffusion map on the PCA
#' library(destiny)
#' dm <- DiffusionMap(t(pca),n_pcs=50, n_eigs = 10)
#' dif.pse <- DPT(dm, tips = c(root, branch.tips), w_width = 0.1)
#' pb.milo$pseudotime <- dif.pse[[DPTroot]]
#' 
#' pb.milo <- markov_probability(
#' milo=pb.milo, 
#' diffusionmap=dm, 
#' terminal_state=branch.tips, 
#' root_cell=root, 
#' pseudotime_key="pseudotime")
#' 
#' #visualization of the results
#' plotPCA(pb.milo,  color_by = "CD8+T") + scale_color_gradientn(colors = pal)
#' plotPCA(pb.milo,  color_by = "CD4+T") + scale_color_gradientn(colors = pal)
#' @export
markov_probability <- function(
    milo, diffusionmap, diffusiontime = NULL, terminal_state, root_cell,
    pseudotime_key = "pseudotime", 
    scale_components = TRUE, num_waypoints = 500) {
  if(is.null(milo[[pseudotime_key]]))
  {
    requireNamespace("rlang")
    if(is.null(diffusiontime))
    {
      rlang::abort(paste("Missing pseudotime data. This data can be either stored in",deparse(substitute(milo)),"or provided by parameter diffusionmap"))
    }
  }
  else
  {
    diffusiontime <- milo[[pseudotime_key]]
  }
  
  # scale data
  multiscale <- .determine.multiscale.space(diffusionmap)
  if (scale_components)
    multiscale <- .minmax.scale(multiscale)
  # sample waypoints to construct markov chain
  waypoints <- .max.min.sampling(multiscale, num_waypoints = 500)
  waypoints <- unique(c(root_cell, waypoints, terminal_state))
  # calculate probabilities
  probabilities <- differentiation_probabilities(multiscale[waypoints, ], terminal_states = terminal_state,
    pseudotime = diffusiontime, waypoints = waypoints)
  # project probabilities from waypoints to each pseudobulk
  probabilities_proj <- project_probability(diffusionmap, waypoints, probabilities)
  # store the result into milo
  requireNamespace("S4Vectors")
  new_coldata <- S4Vectors::DataFrame(probabilities_proj[,
    1], probabilities_proj[, 2])
  colnames(new_coldata) <- c(names(terminal_state))
  # prevent same name in colData
  idx <- names(colData(milo)) %in% colnames(new_coldata)
  if (any(idx)) {
    warning(paste("Name", paste(names(colData(milo))[idx], collapse = ", "),
      "already exists in", as.character(substitute(milo))))
    repeat {
      answer <- readline(prompt = "Do you want to overwrite the column? (y/n): ")
      if (answer == "n") {
        while (any(names(colData(milo)) %in% colnames(new_coldata))) {colnames(new_coldata) <- paste0(colnames(new_coldata), "_new")}
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
