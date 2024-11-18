#' markov_probability
#' 
#' Construct markov chain and calculate probabilities
#' @param milo milo or SingelCellExperiment object, used to store the result
#' @param diffusionmap DiffusionMap object
#' @param scale_components logical
#' If True, the components will be scale
#' @param num_waypoints integer
#' Number of waypoints to sample. 500L by default
#' @param dpt DPT object
#' @param terminal_state 
#' @param root_cell 
#' @return 
#' @export
markov_probability <- function(milo, diffusionmap,dpt,terminal_state, root_cell ,scale_componet = TRUE, num_waypoints = 500  )
{
  
  # scale data 
  multiscale <- .determine.multiscale.space(diffusionmap)
  if(scale_componet) multiscale <- .minmax.scale(multiscale)
  
  # sample waypoints to construct markov chain
  waypoints <- .max.min.sampling(multiscale, num_waypoints=500)
  waypoints <- unique(c(root_cell, waypoints, terminal_state))
  
  # calculate probabilities
  probabilities <- differentiation_probabilities(multiscale[waypoints,],terminal_states = terminal_state, pseudotime = dpt[[paste0("DPT",root_cell)]], waypoints = waypoints) 
  
  # project probabilities from waypoints to each pseudobulk
  probabilities_proj <- project_probability(diffusionmap,waypoints,probabilities)
  
  # store the result into milo
  new_coldata <- DataFrame(dpt[[paste0("DPT",root_cell)]], probabilities_proj[,1], probabilities_proj[,2])
  colnames(new_coldata) <- c("pseudotime", names(terminal_state))
  
  
  # prevent same name in colData
  idx <- names(colData(milo)) %in% colnames(new_coldata)
  if(any(idx))
  {
    warning(paste("Name", paste(names(colData(milo))[idx], collapse = ", "),"already exists in", as.character(substitute(milo))))
    repeat{
      answer <- readline(prompt = "Do you want to overwrite the column? (y/n): ") 
      if(answer == "n")
      {
        while(any(names(colData(milo))%in%colnames(new_coldata))) colnames(new_coldata) <- paste0(colnames(new_coldata),"_new")
        message(paste("The data will stored in", paste(colnames(new_coldata),collapse = ", ")))
        break
      }
      else if(answer == "y")
      {
        message(paste0("Overwriting ", paste(names(colData(milo))[idx], collapse = ", "),"..."))
        colData(milo) <- colData(milo)[!idx]
        break
      }
      else
      {
        message("Invalid response. Please enter 'y' or 'n'.")
      }
    }
  }
  
  
  colData(milo)  <- cbind(colData(milo),new_coldata)
  return(milo)
}
