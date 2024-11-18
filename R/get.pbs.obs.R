#' .get.pbs.obs
#' 
#' Helper function to create the new pseudobulk object's coldata. 
#' @result pbs_obs
.get.pbs.obs <- function(pbs,obs_to_take, milo) 
{
  # prepare per-pseudobulk calls of specified metadata columns
  if(!is.null(obs_to_take))
  {
    pbs.obs <- DataFrame()
    for(anno_col in obs_to_take)
    {
      fa <- as.factor(colData(milo)[[anno_col]])
      fa <- data.frame(fa)
      anno.dummies <- model.matrix(~ . - 1, data = fa, contrasts.arg = lapply(fa, contrasts, contrasts = FALSE))
      anno.count <- t(pbs) %*% anno.dummies   # t(cell x  pseudo)   %*% (cell x element) = pseudo x element
      anno.sum <- apply(anno.count, 1, sum)
      anno.frac<- anno.count/anno.sum
      anno.frac<- DataFrame(anno.frac)
      colnames(anno.frac)<- colnames(anno.dummies)
      pbs.obs[anno_col]<- colnames(anno.frac)[apply(anno.frac,1,which.max)]
      pbs.obs[paste0(anno_col,"_fraction")] <- apply(anno.frac,1,max)
    }
  }
  # report the number of cells for each pseudobulk
  # ensure pbs is an array so that it sums into a vector that can go in easily
  pbs.obs["cell_count"] <- apply(pbs,2, sum)
  return(pbs.obs)
}
