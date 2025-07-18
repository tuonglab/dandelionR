#' Preprocess V(D)J Data for Pseudobulk Analysis
#'
#' This function preprocesses single-cell V(D)J sequencing data for
#' pseudobulk analysis. It filters data based on productivity and chain status,
#' subsets data, extracts main V(D)J genes, and removes unmapped entries.
#'
#' @param sce A `SingleCellExperiment` object. V(D)J data should be contained
#'  in `colData` for filtering.
#' @param mode_option Optional character. Specifies the mode for extracting
#'  V(D)J genes.
#' If `NULL`, `extract_cols` must be specified. Default is `NULL`.
#' @param already.productive Logical. Whether the data has already been filtered
#'  for productivity.
#' If `TRUE`, skips productivity filtering. Default is `FALSE`.
#' @param productive_cols Character vector. Names of `colData` columns used for
#'  productivity filtering.
#' Default is `NULL`.
#' @param productive_vj Logical. If `TRUE`, retains cells where the main
#'  VJ chain is productive.
#'  Default is `TRUE`.
#' @param productive_vdj Logical. If `TRUE`, retains cells where the
#'  main VDJ chain is productive.
#' Default is `TRUE`.
#' @param allowed_chain_status Character vector. Specifies chain statuses to
#'  retain. Valid options
#' include\code{`c('single pair', 'Extra pair', 'Extra pair-exception',
#' 'Orphan VDJ', 'Orphan VDJ-exception')`}. Default is `NULL`.
#' @param subsetby Character. Name of a `colData` column for subsetting.
#'  Default is `NULL`.
#' @param groups Character vector. Specifies the subset condition for filtering.
#'  Default is `NULL`.
#' @param extract_cols Character vector. Names of `colData` columns where V(D)J
#'  information is
#' stored, used instead of the standard columns. Default is `NULL`.
#' @param filter_unmapped Logic. Whether to filter unmapped data. Default
#'  is TRUE.
#' @param check_vj_mapping Logic vector. Whether to check for VJ mapping.
#'  Default is `c(TRUE, TRUE)`.
#'  - If the first element is TRUE, function will filter the unmapped data in V
#'   gene of the VJ chain
#'  - If the second element is TRUE, function will filter the unmapped data in J
#'   gene of the VJ chain
#' @param check_vdj_mapping Logic vector. Specifies columns to check for
#'  VDJ mapping. Default
#' is `c(TRUE, FALSE, 'TRUE)`.
#'  - If the first element is TRUE, function will filter the unmapped data in V
#'   gene of the VDJ chain
#'  - If the second element is TRUE, function will filter the unmapped data in D
#'   gene of the VDJ chain
#'  - If the third element is TRUE, function will filter the unmapped data in J
#'   gene of the VDJ chain
#' @param check_extract_cols_mapping Character vector. Specifies columns related
#'  to `extract_cols`
#' for mapping checks. Default is `NULL`.
#' @param remove_missing Logical. If `TRUE`, removes cells with contigs matching
#'  the filter.
#' If `FALSE`, masks them with uniform values. Default is `TRUE`.
#' @param verbose Logical. Whether to print messages. Default is `TRUE`.

#' @details
#' The function performs the following preprocessing steps:
#' - **Productivity Filtering**:
#'   - Skipped if `already.productive = TRUE`.
#'   - Filters cells based on productivity using `productive_cols` or standard
#'   `colData` columns named `productive_{mode_option}_{type}` (where `type`
#'   is 'VDJ' or 'VJ').
#'   - *mode_option*
#'      - function will check colData(s) named
#'      `productive_{mode_option}_{type}`, where type should be 'VDJ' or 'VJ'
#'      or both, depending on values of productive_vj and productive_vdj.
#'      - If set as `NULl`, the function needs the option 'extract_cols' to be
#'       specified
#'   - *productive_cols*
#'      - must be be specified when productivity filtering is need to conduct
#'       and mode_option is NULL.
#'      - where VDJ/VJ information is stored so that this will be used
#'       instead of the standard columns.
#'   - *productive_vj, productive_vdj*
#'      - If `TRUE`, cell will only be kept if the main V(D)J chain
#'       is productive
#' - **Chain Status Filtering**:
#'   - Retains cells with chain statuses specified by `allowed_chain_status`.
#' - **Subsetting**:
#'   - Conducted only if both `subsetby` and `groups` are provided.
#'   - Retains cells matching the `groups` condition in the `subsetby` column.
#' - **Main V(D)J Extraction**:
#'   - Uses `extract_cols` to specify custom columns for
#'    extracting V(D)J information.
#' - **Unmapped Data Filtering**:
#'   - decided to removes or masks cells based on `filter_unmapped`.
#'   - Checks specific columns for unclear mappings using `check_vj_mapping`,
#'    `check_vdj_mapping`, or `check_extract_cols_mapping`.
#'   - *filter_unmapped*
#'      - pattern to be filtered from object.
#'      - If is set to be `NULL`, the filtering process will not start
#'   - *check_vj_mapping, check_vdj_mapping*
#'      - only `colData` specified by these arguments
#'       (`check_vj_mapping` and `check_vdj_mapping`) will be checked
#'        for unclear mappings
#'   - *check_extract_cols_mapping, related to extract_cols*
#'      - Only `colData` specified by the argument will be checked for
#'      unclear mapping, the colData should first specified by extract_cols
#'   - remove_missing
#'      - If `TRUE`, will remove cells with contigs matching the filter from the
#'       object.
#'      - If `FALSE`, will mask them with a uniform value dependent on
#'      the column name.
#' @include check.R
#' @include filterCells.R
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData<-
#' @importFrom rlang abort
#' @return filtered SingleCellExperiment object
#' @examples
#'
#' # load data
#' data(sce_vdj)
#' # check the dimension
#' dim(sce_vdj)
#' # filtered the data
#' sce_vdj <- setupVdjPseudobulk(
#'     sce = sce_vdj,
#'     mode_option = "abT", # set the mode to alpha-beta TCR
#'     allowed_chain_status = c("Single pair", "Extra pair"),
#'     already.productive = FALSE
#' ) # need to filter the unproductive cells
#' # check the remaining dim
#' dim(sce_vdj)
#'
#' @export
setupVdjPseudobulk <- function(
    sce,
    mode_option = c("abT", "gdT", "B"),
    already.productive = TRUE,
    productive_cols = NULL, productive_vj = TRUE,
    productive_vdj = TRUE,
    allowed_chain_status = NULL, subsetby = NULL,
    groups = NULL, extract_cols = NULL,
    filter_unmapped = TRUE,
    check_vj_mapping = c(TRUE, TRUE),
    check_vdj_mapping = c(TRUE, FALSE, TRUE),
    check_extract_cols_mapping = NULL,
    remove_missing = TRUE,
    verbose = TRUE) {
    # check if the data type is correct
    .classCheck(sce, "SingleCellExperiment")
    mode_option <- match.arg(mode_option)
    .typeCheck(productive_cols, "character")
    .typeCheck(productive_vdj, "logical")
    .typeCheck(productive_vj, "logical")
    .typeCheck(extract_cols, "character")
    .typeCheck(filter_unmapped, "logical")
    ## filter out cells with unproductive chain
    if (!already.productive) {
        sce <- .filterProductivity(
            sce, mode_option, productive_cols,
            productive_vj, productive_vdj, verbose
        )
    }
    ## retain only cells with allowed chain status
    sce <- .allowedChain(sce, allowed_chain_status, verbose)
    ## subset sce by subsetby and groups
    sce <- .subsetSce(sce, subsetby, groups, verbose)
    ## extract main VDJ from specified columns
    if (verbose) message("VDJ data extraction begin:")
    sce <- .extractVdj(sce, extract_cols, mode_option, verbose)
    # remove unclear mapping
    if (filter_unmapped) {
        sce <- .filterUnmapped(
            sce$sce, mode_option, check_vj_mapping,
            check_vdj_mapping,
            sce$main_cols,
            check_extract_cols_mapping, remove_missing, verbose
        )
    } else {
        sce <- sce$sce
    }
    if (verbose) message(sprintf("%d of cells remain.", dim(sce)[2]))
    return(sce)
}

#' filer out cell with unproductive chain
#'
#' @param sce SingleCellExperiment input
#' @param mode_option check setupVdjPseudobulk for detailed explanation
#' @param productive_vj If `TRUE`, retains cells where the main VJ chain is
#'  productive.
#' @param productive_vdj If `TRUE`, retains cells where the main VDJ chain is
#'  productive.
#' @param verbose logical, print messages. Default is `TRUE`.
#' @keywords internal
#' @include check.R
#' @import SingleCellExperiment
#' @importFrom rlang abort
#' @return SingleCellExperiment object after filtering on producive chain
.filterProductivity <- function(sce, mode_option, productive_cols,
                                productive_vj, productive_vdj, verbose) {
    if (is.null(mode_option)) {
        if (!is.null(productive_cols)) { # nocov start
            msg <- paste(productive_cols, collapse = ", ")
            if (verbose) {
                message(sprintf("Checking productivity from %s ..."),
                    appendLF = FALSE
                )
            }
            cnumber0 <- dim(sce)[2]
            sce <- Reduce(function(data, p_col) {
                idx <- substr(colData(data)[[p_col]],
                    start = 1,
                    stop = 1
                ) == "T"
                data[, idx]
            }, productive_cols, init = sce)
            cnumber1 <- dim(sce)[2]
            if (verbose) {
                message(sprintf("%d of cells filtered", cnumber0 - cnumber1))
            }
        } else {
            abort(sprintf("Specify productive_cols %s when mode_option = NULL"))
        } # nocov end
    } else {
        if (is.null(productive_cols)) {
            produ_col <- paste("productive", mode_option, c("VDJ", "VJ"),
                sep = "_"
            )[c(productive_vdj, productive_vj)]
        } else {
            produ_col <- productive_cols
        }
        msg <- paste(produ_col, collapse = ", ")
        if (verbose) {
            message(sprintf("Checking productivity from %s ...", msg),
                appendLF = FALSE
            )
        }
        cnum0 <- dim(sce)[2]
        sce <- Reduce(function(data, p_col) {
            idx <- substr(colData(data)[[p_col]], start = 1, stop = 1) == "T"
            data[, idx]
        }, produ_col, init = sce)
        cnum1 <- dim(sce)[2]
        if (verbose) message(sprintf("%d of cells filtered", cnum0 - cnum1))
    }
    return(sce)
}

#' filtering cell without allowed chain status
#'
#' @param sce SingleCellExperiment object input
#' @param allowed_chain_status the chain needs to be retain, passed from
#'  setupVdjPseudobulk
#' @param verbose logical, print messages. Default is `TRUE`.
#' @keywords internal
#' @import SingleCellExperiment
#' @importFrom rlang abort
#' @return SingleCellExperiment object with allowed chain status
.allowedChain <- function(sce, allowed_chain_status, verbose) {
    ## retain only cells with allowed chain status
    if (!is.null(allowed_chain_status)) {
        if (verbose) message("checking allowed chains...", appendLF = FALSE)
        cnumber0 <- dim(sce)[2]
        idx <- colData(sce)[["chain_status"]] %in% allowed_chain_status
        if (!any(idx)) {
            allowed_cs <- paste(allowed_chain_status, collapse = ", ")
            current_cs <- paste(unique(colData(sce)[["chain_status"]]),
                collapse = ", "
            )
            abort(sprintf(
                paste0(
                    "Unsuitable allowed_chain_status,\n ",
                    "The current allowed_chain_status: %s.\n ",
                    "While the chain status in the dataset: %s."
                ),
                allowed_cs, current_cs
            ))
        }
        sce <- sce[, idx]
        cnumber1 <- dim(sce)[2]
        filtered <- cnumber0 - cnumber1
        if (verbose) message(sprintf("%d of cells filtered", filtered))
    }
    return(sce)
}
#' Subset sce with given parameter
#'
#' @param sce SingleCellExperiment object input
#' @param subsetby subsetby Character. Name of a `colData` column for
#'  subsetting. given by setupVdjPsudobulk.
#' @param groups Character vector. Specifies the subset condition for filtering.
#'  given by setupVdjPsudobulk.
#' @param verbose logical, print messages. Default is `TRUE`.
#' @keywords internal
#' @import SingleCellExperiment
#' @importFrom rlang abort
#' @return subsetted SingleCellExperiment object
.subsetSce <- function(sce, subsetby, groups, verbose) {
    .typeCheck(subsetby, "character")
    .typeCheck(groups, "character")
    if (!is.null(groups) && !is.null(subsetby)) {
        msg1 <- paste(as.character(substitute(groups))[-1], collapse = ", ")
        msg2 <- as.character(substitute(subsetby))
        if (verbose) {
            message(sprintf("Subsetting data with %s in %s ...", msg1, msg2),
                appendLF = FALSE
            )
        }
        cnumber0 <- dim(sce)[2]
        idx <- Reduce(`|`, lapply(groups, function(i) {
            colData(sce)[[subsetby]] %in% i
        }))
        sce <- sce[, idx]
        cnumber1 <- dim(sce)[2]
        filtered <- cnumber0 - cnumber1
        if (verbose) message(sprintf("%d of cells filtered", filtered))
    }
    return(sce)
}
#' Specify the columns which store VDJ information, and extract the main chain
#'  from it
#'
#' @param sce SingleCellExperiment object input
#' @param extract_cols The setupVdjPseutobulk transfered parameter given by user
#'  to specify the VDJ information columns
#' @param mode_option see document of setupVdjPseudobulk for
#'  detailed explanation
#' @param verbose logical, print messages. Default is `TRUE`.
#' @keywords internal
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData<-
#' @importFrom rlang abort
#' @return SingleCellExperiment objects with column stores the information of
#'  the main VDJ information in colData slot
.extractVdj <- function(sce, extract_cols, mode_option, verbose) {
    # generate colnames to extract
    if (is.null(extract_cols)) {
        extract_cols <- .generateExtractName(sce, mode_option, verbose)
    }
    # make sure the columns exist
    sce <- .generateExtractColumn(sce, extract_cols, verbose)
    # extract main information
    if (!length(grep("_VDJ_main|_VJ_main", names(colData(sce))))) {
        colns <- paste(extract_cols, collapse = ", ")
        if (verbose) {
            message(sprintf("Extract main TCR from %s ...", colns),
                appendLF = FALSE
            )
        }
        sce <- Reduce(function(data, ex_col) {
            tem <- colData(data)[[ex_col]]
            strtem <- strsplit(as.character(tem), "\\|")
            value_need <- vapply(strtem, `[`, 1, FUN.VALUE = character(1))
            colData(data)[[paste(ex_col, "main", sep = "_")]] <- value_need
            data
        }, extract_cols, init = sce)
        main_cols <- paste(extract_cols, "main", sep = "_")
        if (verbose) message("Complete.")
    } else {
        main_cols <- colnames(colData(sce))[grep(
            "_VDJ_main|_VJ_main",
            names(colData(sce))
        )]
        if (verbose) {
            warning(
                c(
                    "main VDJ information already exists,",
                    "Instead of using 'check_v(d)j_mapping', please use",
                    "the argument 'check_extract_cols_mapping'",
                    "to clarify the columns undergo filtering."
                )
            )
        }
    }
    return(list(sce = sce, main_cols = main_cols))
}
#' Generate the name of columns with given parameter
#'
#' @param sce SingleCellExperiment object input
#' @param mode_option see document of setupVdjPseudobulk for explanation
#' @param verbose logical, print messages. Default is `TRUE`.
#' @keywords internal
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData<-
#' @return a vecotor of colnames we need to perform main chain extraction
.generateExtractName <- function(sce, mode_option, verbose) {
    if (verbose) {
        message(c(
            "extract_cols not specified, automatically generate",
            " colnames for extraction."
        ))
    }
    v_call <- if ("v_call_genotyped_VDJ" %in% colnames(colData(sce))) {
        # nocov start
        "v_call_genotyped_" # nocov end
    } else {
        "v_call_"
    }
    prefix <- c(v_call, "d_call_", "j_call_")
    if (!is.null(mode_option)) {
        # can be pack as another function
        suffix <- c("_VJ", "_VDJ")
        extr_cols <- as.vector(outer(prefix, suffix, function(x, y) {
            paste0(x, mode_option, y)
        }))
        extr_cols <- extr_cols[extr_cols != paste0(
            "d_call_", mode_option,
            "_VJ"
        )]
    } else {
        # nocov start
        suffix <- c("VJ", "VDJ")
        extr_cols <- as.vector(outer(prefix, suffix, function(x, y) {
            paste0(x, y)
        }))
        extr_cols <- extr_cols[extr_cols != paste0("d_call_", "VJ")]
        # nocov end
    }
    return(extr_cols)
}
#' Check whether the columns with specified names exist, if not, create them
#' with CTgene columns
#' @param sce SingleCellExperiment object input
#' @param extract_cols column names we aim to extract information from
#' @param verbose logical, print messages. Default is `TRUE`.
#' @keywords internal
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData<-
#' @importFrom rlang abort
#' @return SingleCellExperiment with columns containing VDJ information in the
#' names we've specified.
.generateExtractColumn <- function(sce, extract_cols, verbose) {
    msg <- paste(extract_cols, collapse = ", ")
    if (!any(extract_cols %in% colnames(colData(sce)))) {
        if (verbose) {
            message(sprintf(
                "ColData does not exist, Creating %s %s",
                msg, "colData based on column CTgene"
            ))
        }
        if (!"CTgene" %in% colnames(colData(sce))) {
            abort(sprintf(paste(
                "Both %s and CTgene do not exist\n You could modify parameter",
                "extract_cols to clarify VDJ information's location"
            ), msg))
        }
        splitVdj <- splitCTgene(sce)
        if (length(splitVdj[[1]]) != length(extract_cols)) {
            abort(
                sprintf(paste(
                    "Colnames %s must have the same length with the",
                    "vdj data, which is of length %d.\n Modify",
                    "`extract_cols`` to specify the TCR columns."
                )),
                paste0(extract_cols[!extract_cols %in% colnames(colData(sce))],
                    collapse = ", "
                ),
                length(splitVdj[[1]])
            )
        } else {
            vdj <- lapply(seq(length(extract_cols)), function(X, sc) {
                vapply(X = sc, "[", X, FUN.VALUE = character(1))
            }, sc = splitVdj)
            names(vdj) <- extract_cols
            colData(sce) <- cbind(colData(sce), vdj)
        }
    } else if (!all(extract_cols %in% colnames(colData(sce)))) {
        abort(sprintf(
            paste0(
                "Colnames %s do not exist in colData.\n Use",
                "`extract_cols` to specify the TCR columns."
            ),
            paste0(
                extract_cols[!extract_cols %in% colnames(colData(sce))],
                collapse = ", "
            )
        ))
    }
    return(sce)
}
#' Filter out cell with unclear mapping in VDJ information
#' @param sce SingleCellExperiment object input
#' @param mode_option see document of setupVdjPseudobulk for explanation
#' @param check_vj_mapping logical vector to set whether to check V and J gene
#' in VJ chain, passed from setupVdjPseudobulk
#' @param check_vdj_mapping logical vector to set whether to check
#' V, D and J gene in VDJ chain, passed from setupVdjPseudobulk
#' @param main_cols column names in colData in which The information of
#' main chain stores
#' @param check_extract_cols_mapping character vector,the names of columns that
#'  needs to be checked, passed from setupVdjPseudobulk
#' @param remove_missing option for removing the unclear mappin or just mask it,
#' passed from setupVdjPseudobulk
#' @param verbose logical, print messages.
#' @keywords internal
#' @include filterCells.R
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData<-
#' @importFrom rlang abort
#' @return filtered SingleCellExperiment object
.filterUnmapped <- function(sce, mode_option, check_vj_mapping,
                            check_vdj_mapping, main_cols,
                            check_extract_cols_mapping, remove_missing,
                            verbose) {
    .typeCheck(check_vj_mapping, "logical")
    if (length(check_vj_mapping) != 2) {
        abort(sprintf(
            paste0(
                "ValueError: length of check_vj_mapping should be",
                " 2. But %d was provided."
            ),
            length(check_vj_mapping)
        )) # nocov
    }
    .typeCheck(check_vdj_mapping, "logical")
    if (length(check_vdj_mapping) != 3) {
        abort(sprintf(paste0(
            "ValueError: length of check_vj_mapping should be",
            " 3. But %d was provided."
        ), length(check_vdj_mapping))) # nocov
    }
    .typeCheck(check_extract_cols_mapping, "character")
    .typeCheck(remove_missing, "logical")
    filter_pattern <- ",|None|No_contig"
    extr_cols <- c()
    if (is.null(check_extract_cols_mapping)) {
        extr_cols <- main_cols[c(check_vdj_mapping, check_vj_mapping)]
    } else {
        extr_cols <- check_extract_cols_mapping
    }
    if (!is.null(extr_cols)) {
        msg <- paste(extr_cols, collapse = ", ")
        if (verbose) {
            message(sprintf("Filtering cells from %s ...", msg),
                appendLF = FALSE
            )
        }
        cnumber0 <- dim(sce)[2]
        sce <- Reduce(function(x, y) {
            .filterCells(
                sce = x, col_n = y, filter_pattern = filter_pattern,
                remove_missing = remove_missing
            )
        }, extr_cols, init = sce)
        cnumber1 <- dim(sce)[2]
        filtered <- cnumber0 - cnumber1
        if (verbose) message(sprintf("%d of cells filtered", filtered))
    }
    return(sce)
}
