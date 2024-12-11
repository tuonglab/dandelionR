#' Split the V(D)J genes from `CTgene` column and store them separately.
#'
#' @param sce SingleCellExperiment object after conducting scRepertoire::combineTCR()
#' @keywords internal
#' @return list contain vector of VJ & VDJ genes from each cell
#' @importFrom SummarizedExperiment colData<-
splitCTgene <- function(sce) {
    CTgene <- colData(sce)[, "CTgene"]
    # split with '_'
    CTgene <- lapply(CTgene, strsplit, split = "_")
    # split with '.'
    CTgene <- lapply(CTgene, function(x) lapply(x, strsplit, split = "[.]"))
    CTgene <- .collapse_nested_list(CTgene)
    # names(CTgene) <- NULL # this is a list containing `ncol()` sublists,
    # where each sublist contains two character vectors: the first represent
    # TRA the second represent TRB
    CTgene <- lapply(CTgene, formatVdj)
    return(lapply(CTgene, unlist))
}

#' Change the format of splitCTgene output.
#'
#' @param gene_list list containing the output from splitCTgene.
#' @keywords internal
#' @return list contain vector of VJ + VDJ information of the cell input
formatVdj <- function(gene_list) {
    if (all(is.na(gene_list))) {
        vj <- rep("None", 2)
        vdj <- rep("None", 3)
    } else {
        vj <- chainAssign(gene_list[[1]], 2)
        vdj <- chainAssign(gene_list[[2]], 3)
    }
    return(append(vj, vdj))
}


#' Assign the V(D)J gene to the right chain.
#'
#' @param vec vector of V(D)J genes to assign to the right chain.
#' @param num number of genes to return. should be 2(vj) or 3(vdj)
#' @keywords internal
#' @return list contain vector of VJ + VDJ of the cell input
#' @importFrom rlang abort
chainAssign <- function(vec, num) {
    if (length(vec) == (num + 1)) {
        chains <- vec[seq_len(num)]
    } else if (all(vec == "NA")) {
        chains <- rep("None", num)
    } else {
        abort
    }
    return(chains)
}

#' Collapse a nested list
#'
#' @param input_list input nested list.
#' @keywords internal
#' @return collapsed list
.collapse_nested_list <- function(input_list) {
    lapply(input_list, function(sublist) {
        if (is.list(sublist)) {
            # Check if all elements are also lists
            all_are_lists <- TRUE
            for (item in sublist) {
                if (!is.list(item)) {
                    all_are_lists <- FALSE
                    break
                }
            }
            # If all elements are lists, unlist one level
            if (all_are_lists) {
                return(unlist(sublist, recursive = FALSE))
            }
        }
        return(sublist)
    })
}
