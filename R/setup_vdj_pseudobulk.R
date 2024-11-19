#' setup_vdj_pseudobulk
#' 
#' filtering the data based on productivity, chain status, and subset data. ....
#' @param sce SingleCellExperiment object, 
#' vdj data should contained in colData for filtering
#' @param already.productive logical. Whether the data is already filtered by productivity. 
#'  - If true, the function will skip the productivity filtering
#'  - If false, retain only productive entries in columns specified by mode_option or productive_cols
#' @param mode_option optional
#' mode for extraction the V(D)J genes.
#'  - for productive filtering, function will check colData(s) named productive_{mode_option}_{type}, where type should be "VDJ" or "VJ" or both, depending on productive_vj and productive_vdj.
#'    - If set as "NULL", the function needs the option "extract_cols" to be specified
#'  - for 
#' @param productive_cols paramter used in productivity filtering, NULL by default
#'  - colData names where VDJ/VJ information is stored so that this will be used instead of the standard columns.
#'  - must be be specified when productivity filtering is need to conduct and mode_option is NULL. 
#' @param productive_vdj logical, TRUE by default
#' Option in productivity filtering. If True, cell will only be kept if the main VDJ chain is productive
#' @param productive_vj logical, TRUE by default
#' Option in productivity filtering. If True, cell will only be kept if the main VJ chain is productive
#' @param allowed_chain_status character vectors, optional
#' element should come from c("single pair","Extra pair", "Extra pair-exceptipn", "Orphan VDJ","Orphan VDJ-exception")
#' @param subsetby character, NULL by default
#'  - name of one colData provided for subsetting. The subsetting process will only start when both groups and subsetby are provided.
#'  - If provided, only the cell with {subsetby} feature in {groups} will be used for computing the VDJ feature space
#' @param groups character vector, NULL by default
#'  - condition for subsetting
#'  - If provided, only the following groups/ categories will be used for computing the VDJ feature space
#' @param extract_cols character vector, NULL by default
#' column names where VDJ/VJ information is stored so that this will be used instead of the standard columns
#' @param filter_pattern character string, optional ",|None|No_contig" by default
#' pattern to filter form object. If NULL, does not filter
#' @param check_vdj_mapping character vector, optional
#' elements should come from c("v_call", "j_call")
#' Only columns in the argument will be checked for unclear mapping(contain comma) in VDJ calls
#' @param check_vj_mapping character vector, optional
#' elements should come from c("v_call", "j_call")
#' Only columns in the argument will be checked for unclear mapping(contain comma) in VJ calls 
#' @param check_extract_cols_mapping character vecter, NULL by default
#' Only columns in the argument will be checked for unclear mapping (containing comma) in columns specified in extract_cols
#' @param remove_missing bool, True by default
#' If true, will remove cells with contigs matching the filter from the object. If False, will mask them with a uniform value dependent on the column name.
#' @include check.R
#' 
#' @return filtered SingleCellExperiment object
#' @export
setup_vdj_pseudobulk<-function(sce,
         mode_option = c("abT","gdT","B"),
         productive_cols = NULL,
         productive_vdj = TRUE,
         productive_vj = TRUE,
         subsetby = NULL,
         groups = NULL,
         allowed_chain_status = c("Single pair", "Extra pair", "Extra pair-exception", "Orphan VDJ", "Orphan VDJ-exception"),
         extract_cols = NULL,
         filter_pattern = ",|None|No_cotig",
         check_vdj_mapping = c("v_call","j_call"),
         check_vj_mapping = c("v_call","j_call"),
         check_extract_cols_mapping = NULL,
         remove_missing = TRUE,
         already.productive = TRUE,
         ...
)
{
  # check if the data type is correct
  .type.check(sce,"SingleCellExperiment")
  mode_option <- match.arg(mode_option)
  .type.check(productive_cols,"character")
  .type.check(productive_vdj,"logical")
  .type.check(productive_vj,"logical")
  .type.check(subsetby,"character")
  .type.check(groups,"character")
  allowed_chain_status <- match.arg(allowed_chain_status,several.ok = T)
  .type.check(extract_cols,"character")
  .type.check(filter_pattern,"character")
  check_vdj_mapping <- match.arg(check_vdj_mapping,c("v_call","d_call","j_call"), several.ok = T)
  check_vj_mapping <- match.arg(check_vj_mapping, several.ok = T)
  .type.check(check_extract_cols_mapping,"character")
  .type.check(remove_missing,"logical")
  
  
  
  
  # filtering
  ## retain only productive entries based on specified mode
  message("Productivity filtering...")
  if(!already.productive)
  {
    if(is.null(mode_option))
    {
      if(!is.null(productive_cols))
      {
        message(paste("checking column(s) from", paste(productive_cols, collapse = ", ")))
        sce <- Reduce(function(data, p_col){
          idx<- substr(colData(data)[[p_col]],start = 1, stop = 1) == "T"
          data[,idx]
        }, productive_cols, init = sce)
      }
      else
      {
        stop("When mode_option is NULL, the productive_cols must be specified.")
      }
    }
    else
    {
      produ_col <- paste("productive",mode_option,c("VDJ","VJ"),sep = "_")[c(productive_vdj,productive_vj)]
      message(paste("checking column(s) from", paste(produ_col, collapse = ", ")))
      sce <- Reduce(function(data, p_col){
        idx <- substr(colData(data)[[p_col]],start = 1, stop = 1) == "T"
        data[,idx]
      },produ_col, init = sce)
    }
  }

  
  ## retain only cells with allowed chain status
  if(!is.null(allowed_chain_status))
  {
    message("checking allowed chain status")
    idx <- colData(sce)[["chain_status"]] %in% allowed_chain_status
    if(!any(idx))
    {
      abort(paste("Unsuitable allowed_chain_status,\n The current allowed_chain_status:", paste(allowed_chain_status,collapse = ", "), ".\n While the chain status in the dataset:", paste(unique(colData(sce)[["chain_status"]]),collapse = ", ")))
    }
    sce <- sce[,idx]
  }
  
  
  ## subset sce by subsetby and groups
  if (!is.null(groups)&&!is.null(subsetby)){
    idx <- Reduce(`|`, lapply(groups, function(i) colData(sce)[[subsetby]] %in% i))
    sce <- sce[, idx]
  }  
  
  
  ## extract main VDJ from specified columns
  if (is.null(extract_cols))
  {
    if(!length(grep("_VDJ_main|_VJ_main",names(colData(sce))))) 
    {
      v_call<- if("v_call_genotyped_VDJ"%in% colnames(colData(sce))) "v_call_genotyped_" else "v_call_"
      prefix <- c(v_call, "d_call_", "j_call_")
      if(!is.null(mode_option)) # can be pack as another function
      {
        suffix <- c("_VDJ", "_VJ")
        extr_cols <- as.vector(outer(prefix, suffix, function(x,y) paste0(x, mode_option, y)))
        extr_cols <- extr_cols[extr_cols != paste0("d_call_",mode_option, "_VJ")]
        sce <- Reduce(function(data, ex_col){
          tem <- colData(data)[[ex_col]]
          strtem <- strsplit(as.character(tem),"\\|")
          colData(data)[[paste(ex_col,"main",sep = "_")]]<-vapply(strtem, `[`,1, FUN.VALUE = character(1))
          data
        }, extr_cols,init = sce)
        
      }
      else{
        suffix <- c("VDJ", "VJ")
        extr_cols <- as.vector(outer(prefix, suffix, function(x,y) paste0(x,y)))
        extr_cols <- extract_cols[extr_cols != paste0("d_call_", "VJ")]
        sce <- Reduce(function(data, ex_col){
          tem <- colData(data)[[ex_col]]
          strtem <- strsplit(as.character(tem),"\\|")
          colData(data)[[paste(ex_col,"main",sep = "_")]]<-vapply(strtem, `[`,1, FUN.VALUE = character(1))
          data
        }, extr_cols,init = sce)        
      }
    }
  }
  else
  {
    sce <- Reduce(function(data, ex_col){
      tem <- colData(data)[[ex_col]]
      strtem <- strsplit(as.character(tem),"\\|")
      colData(data)[[paste(ex_col,"main",sep = "_")]]<-vapply(strtem, `[`,1, FUN.VALUE = character(1))
      data
    }, extract_cols,init = sce)     
  }
  
  
  # remove unclear mapping
  if(!is.null(filter_pattern))
  {
    extr_cols <- c()
    if(!is.null(mode_option))
    {
      if(!is.null(check_vdj_mapping)) {
        extr_cols <- c(extr_cols, paste(check_vdj_mapping, mode_option,"VDJ_main",sep = "_"))
      }
      if(!is.null(check_vj_mapping)) {
        extr_cols <- c(extr_cols, paste(check_vj_mapping, mode_option,"VJ_main",sep = "_"))
      }
    }
    else
    {
      if(!is.null(extract_cols)) 
      {
        if(!is.null(check_vdj_mapping)) {
          extr_cols <- c(extr_cols, paste(check_vdj_mapping,"VDJ_main",sep = "_"))
        }
        if(!is.null(check_vj_mapping)) {
          extr_cols <- c(extr_cols, paste(check_vj_mapping,"VJ_main",sep = "_"))}
      }
      else
      {
        if(is.null(check_extract_cols_mapping)) extr_cols <- check_extract_cols_mapping
      }
    }
    if(!is.null(extr_cols))
      sce <- Reduce(function(x, y) {
        .filter.cells(sce = x, col_n = y, filter_pattern = filter_pattern, remove_missing = remove_missing)
      },
      extr_cols,
      init = sce
      )
  }
  return(sce)  
  
  
}

