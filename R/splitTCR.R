#' split_CTgene
#'
#' split the vdj gene from CTgene and store them separately in the columns
#' @param combined_Expression SingleCellExperiment object after conducting scRepertoire::combineTCR()
#' @keywords internal
#' @return list contain vector of vj & vdj of TCR from each cell
#' @import SingleCellExperiment
splitCTgene <- function(combined_Expression) {
    CTgene <- colData(combined_Expression)["CTgene"]
    # split with "_"
    CTgene <- sapply(CTgene, strsplit, split = "_")
    # split with "."
    CTgene <- sapply(CTgene, strsplit, split = "[.]")
    names(CTgene) <- NULL # this is a list containing `ncol()` sublists, where each sublist contains two character vectors: the first represent TRA the second represent TRB
    CTgene <- lapply(CTgene, categoryTCR)
    return(lapply(CTgene, unlist))
}

#' categoryTCR
#'
#' change the format of splited CTgene
#' @param gene_list list which is the result of CTgene split
#' @keywords internal
#' @return list contain vector of TCRA + TCRB of the cell input
#' @import SingleCellExperiment
categoryTCR <- function(gene_list) {
    if (all(is.na(gene_list))) {
        TCRA <- rep("None", 2)
        TCRB <- rep("None", 3)
    } else {
        TCRA <- vecAssign(gene_list[[1]], 2)
        TCRB <- vecAssign(gene_list[[2]], 3)
    }
    return(append(TCRA, TCRB))
}


#' vecAssign
#'
#' Assgin the vdj(vj) gene to the result
#' @param vec list of vdj(vj) gene to assign to the reusult
#' @param num number of gene should return, should be 2(vj) or 3(vdj)
#' @keywords internal
#' @return list contain vector of TCRA + TCRB of the cell input
#' @import SingleCellExperiment
#' @importFrom rlang abort
vecAssign <- function(vec, num) {
    if (length(vec) == (num + 1)) {
        TCR <- vec[1:num]
    } else if (all(vec == "NA")) {
        TCR <- rep("None", num)
    } else {
        abort
    }
    return(TCR)
}
