#'.get.pbs
#' 
#'Helper function to ensure we have cells by pseudobulks matrix which we can use for pseudobulking. Use the pbs and obs_to_bulk inputs to vdj_pseudobulk() and gex_pseudobulk()
#'
.get.pbs<-function(pbs, obs_to_bulk, milo)
{
  # some way to pseudobulk
  if(is.null(pbs)&&is.null(obs_to_bulk))
  {
    abort("You need to specify 'pbs' or 'obs_to_bulk', not both")
  }
  # but just one
  if(!is.null(pbs)&&!is.null(obs_to_bulk))
  {
    abort("You need to specify 'pbs' or 'obs_to_bulk', not both" )
  }
  # turn the pesudobulk matrix dense if need be?
  if(!is.null(pbs))
  {
    return(pbs)
  }
  if(!is.null(obs_to_bulk))
  {
    tobulk <- lapply(obs_to_bulk,function(x){
      colData(milo)[[x]]
    })
    names(tobulk) <- obs_to_bulk
    tobulk <- as.data.frame(tobulk)
    tobulk<-as.data.frame(apply(tobulk, 1, paste,collapse = ",",simplify = FALSE))
    tobulk <- model.matrix(~ t(tobulk) - 1)
    colnames(tobulk) <- gsub("t\\(tobulk\\)","",colnames(tobulk))
    return(tobulk)
  }
}
