#' markov_probability
#'
#' Preprocessing data and Construct markov chain and calculate probabilities
#' @param milo milo or SingelCellExperiment object, used to store the result
#' @param diffusionmap DiffusionMap object corresponds to milo
#' @param scale_components logical, If True, the components will be scale
#' @param num_waypoints integer, 500L by default. 
#' - Number of waypoints to sample to construct markov chain. 
#' @param diffusiontime numeric vector, contain the difussion time of each pseudobulk
#' @param terminal_state the index of the terminal state
#' @param root_cell the index of the root state
#' @return milo or SinglCellExperiment object with pseudotime, probabilities in its colData
#' @details
#' The preprocessing process contain data scaling and waypoints selecting
#' 
#' @include determ.multiscale.space.R
#' @export
markov_probability <- function(milo, diffusionmap, diffusiontime, terminal_state, root_cell,
  scale_componet = TRUE, num_waypoints = 500) {
  # scale data
  multiscale <- .determine.multiscale.space(diffusionmap)
  if (scale_componet)
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
  new_coldata <- DataFrame(diffusiontime, probabilities_proj[,
    1], probabilities_proj[, 2])
  colnames(new_coldata) <- c("pseudotime", names(terminal_state))


  # prevent same name in colData
  idx <- names(colData(milo)) %in% colnames(new_coldata)
  if (any(idx)) {
    warning(paste("Name", paste(names(colData(milo))[idx], collapse = ", "),
      "already exists in", as.character(substitute(milo))))
    repeat {
      answer <- readline(prompt = "Do you want to overwrite the column? (y/n): ")
      if (answer == "n") {
        while (any(names(colData(milo)) %in% colnames(new_coldata))) colnames(new_coldata) <- paste0(colnames(new_coldata),
          "_new")
        message(paste("The data will stored in", paste(colnames(new_coldata),
          collapse = ", ")))
        break
      } else if (answer == "y") {
        message(paste0("Overwriting ", paste(names(colData(milo))[idx], collapse = ", "),
          "..."))
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