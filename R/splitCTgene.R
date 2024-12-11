#' Split the V(D)J genes from `CTgene` column and store them separately.
#'
#' @param sce SingleCellExperiment object after conducting scRepertoire::combineTCR()
#' @keywords internal
#' @return list contain vector of VJ & VDJ genes from each cell
#' @importFrom SummarizedExperiment colData<-
splitCTgene <- function(sce) {
    CTgene <- colData(sce)["CTgene"]
    # split with "_"
    CTgene <- vapply(CTgene, strsplit, split = "_", FUN.VALUE = list(character(1)))
    # split with "."
    CTgene <- vapply(CTgene, strsplit, split = "[.]", FUN.VALUE = list(character(1)))
    names(CTgene) <- NULL # this is a list containing `ncol()` sublists, where each sublist contains two character vectors: the first represent TRA the second represent TRB
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
