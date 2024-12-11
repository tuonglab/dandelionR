#' Preprocess V(D)J Data for Pseudobulk Analysis
#'
#' This function preprocesses single-cell V(D)J sequencing data for pseudobulk analysis. It filters data based
#' on productivity and chain status, subsets data, extracts main V(D)J genes, and removes unmapped entries.
#'
#' @param sce A `SingleCellExperiment` object. V(D)J data should be contained in `colData` for filtering.
#' @param mode_option Optional character. Specifies the mode for extracting V(D)J genes.
#' If `NULL`, `extract_cols` must be specified. Default is `NULL`.
#' @param already.productive Logical. Whether the data has already been filtered for productivity.
#' If `TRUE`, skips productivity filtering. Default is `FALSE`.
#' @param productive_cols Character vector. Names of `colData` columns used for productivity filtering.
#' Default is `NULL`.
#' @param productive_vj Logical. If `TRUE`, retains cells where the main VJ chain is productive.
#'  Default is `TRUE`.
#' @param productive_vdj Logical. If `TRUE`, retains cells where the main VDJ chain is productive.
#' Default is `TRUE`.
#' @param allowed_chain_status Character vector. Specifies chain statuses to retain. Valid options
#' include `c('single pair', 'Extra pair', 'Extra pair-exception', 'Orphan VDJ', 'Orphan VDJ-exception')`. Default is `NULL`.
#' @param subsetby Character. Name of a `colData` column for subsetting. Default is `NULL`.
#' @param groups Character vector. Specifies the subset condition for filtering. Default is `NULL`.
#' @param extract_cols Character vector. Names of `colData` columns where V(D)J information is
#' stored, used instead of the standard columns. Default is `NULL`.
#' @param filter_unmapped Logic. Whether to filter unmapped data. Default is TRUE.
#' @param check_vj_mapping Logic vector. Whether to check for VJ mapping. Default
#' is `c(TRUE, TRUE)`. 
#'  - If the first element is TRUE, function will filter the unmapped data in V gene of the VJ chain
#'  - If the second element is TRUE, function will filter the unmapped data in J gene of the VJ chain
#' @param check_vdj_mapping Logic vector. Specifies columns to check for VDJ mapping. Default
#' is `c(TRUE, FALSE, 'TRUE)`.
#'  - If the first element is TRUE, function will filter the unmapped data in V gene of the VDJ chain 
#'  - If the second element is TRUE, function will filter the unmapped data in D gene of the VDJ chain   
#'  - If the third element is TRUE, function will filter the unmapped data in J gene of the VDJ chain   
#' @param check_extract_cols_mapping Character vector. Specifies columns related to `extract_cols`
#' for mapping checks. Default is `NULL`.
#' @param remove_missing Logical. If `TRUE`, removes cells with contigs matching the filter.
#' If `FALSE`, masks them with uniform values. Default is `TRUE`.
#'
#' @details
#' The function performs the following preprocessing steps:
#' - **Productivity Filtering**:
#'   - Skipped if `already.productive = TRUE`.
#'   - Filters cells based on productivity using `productive_cols` or standard `colData` columns named `productive_{mode_option}_{type}` (where `type` is 'VDJ' or 'VJ').
#'   - *mode_option*
#'      - function will check colData(s) named `productive_{mode_option}_{type}`, where type should be 'VDJ' or 'VJ' or both, depending on values of productive_vj and productive_vdj.
#'      - If set as `NULl`, the function needs the option 'extract_cols' to be specified
#'   - *productive_cols*
#'      - must be be specified when productivity filtering is need to conduct and mode_option is NULL.
#'      - where VDJ/VJ information is stored so that this will be used instead of the standard columns.
#'   - *productive_vj, productive_vdj*
#'      - If `TRUE`, cell will only be kept if the main V(D)J chain is productive
#' - **Chain Status Filtering**:
#'   - Retains cells with chain statuses specified by `allowed_chain_status`.
#' - **Subsetting**:
#'   - Conducted only if both `subsetby` and `groups` are provided.
#'   - Retains cells matching the `groups` condition in the `subsetby` column.
#' - **Main V(D)J Extraction**:
#'   - Uses `extract_cols` to specify custom columns for extracting V(D)J information.
#' - **Unmapped Data Filtering**:
#'   - decided to removes or masks cells based on `filter_unmapped`.
#'   - Checks specific columns for unclear mappings using `check_vj_mapping`, `check_vdj_mapping`, or `check_extract_cols_mapping`.
#'   - *filter_unmapped*
#'      - pattern to be filtered from object.
#'      - If is set to be `NULL`, the filtering process will not start
#'   - *check_vj_mapping, check_vdj_mapping*
#'      - only `colData` specified by these arguments (`check_vj_mapping` and `check_vdj_mapping`) will be checked for unclear mappings
#'   - *check_extract_cols_mapping, related to extract_cols*
#'      - Only `colData` specified by the argument will be checked for unclear mapping, the colData should first specified by extract_cols
#'   - remove_missing
#'      - If `TRUE`, will remove cells with contigs matching the filter from the object.
#'      - If `FALSE`, will mask them with a uniform value dependent on the column name.
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
#' allowed <- c("Single pair", "Extra pair", "Extra pair-exception", "Orphan VDJ", "Orphan VDJ-exception")
#' # filtered the data
#' sce_vdj <- setupVdjPseudobulk(
#'     sce = sce_vdj,
#'     mode_option = "abT", # set the mode to αβTCR
#'     allowed_chain_status = allowed,
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
    productive_cols = NULL, productive_vj = TRUE, productive_vdj = TRUE, allowed_chain_status = NULL,
    subsetby = NULL, groups = NULL, extract_cols = NULL, filter_unmapped = TRUE,
    check_vj_mapping = c(TRUE, TRUE), check_vdj_mapping = c(TRUE, FALSE, TRUE),
    check_extract_cols_mapping = NULL, remove_missing = TRUE) {
    # check if the data type is correct
    .classCheck(sce, "SingleCellExperiment")
    mode_option <- match.arg(mode_option)
    .typeCheck(productive_cols, "character")
    .typeCheck(productive_vdj, "logical")
    .typeCheck(productive_vj, "logical")
    .typeCheck(subsetby, "character")
    .typeCheck(groups, "character")
    .typeCheck(extract_cols, "character")
    .typeCheck(filter_unmapped, "logical")
    .typeCheck(check_vj_mapping, "logical")
    if(length(check_vj_mapping)!=2)
        abort(paste("ValueError: length of check_vj_mapping should be 2. But",length(check_vj_mapping),"was provided."))
    .typeCheck(check_vdj_mapping, "logical")
    if(length(check_vdj_mapping)!=3)
        abort(paste("ValueError: length of check_vj_mapping should be 3. But",length(check_vdj_mapping),"was provided."))
    .typeCheck(check_extract_cols_mapping, "character")
    .typeCheck(remove_missing, "logical")

    # filtering retain only productive entries based on specified mode
    if (!already.productive) {
        if (is.null(mode_option)) {
            if (!is.null(productive_cols)) {
                msg <- paste(productive_cols, collapse = ", ")
                message(sprintf("Checking productivity from %s ..."), appendLF = FALSE)
                cnumber0 <- dim(sce)[2]
                sce <- Reduce(function(data, p_col) {
                    idx <- substr(colData(data)[[p_col]], start = 1, stop = 1) == "T"
                    data[, idx]
                }, productive_cols, init = sce)
                cnumber1 <- dim(sce)[2]
                filtered <- cnumber0 - cnumber1
                message(sprintf("%d of cells filtered", filtered))
            } else {
                abort("When mode_option is NULL, the productive_cols must be specified.")
            }
        } else {
            produ_col <- paste("productive", mode_option, c("VDJ", "VJ"), sep = "_")[c(
                productive_vdj,
                productive_vj
            )]
            msg <- paste(produ_col, collapse = ", ")
            message(sprintf("Checking productivity from %s ...", msg), appendLF = FALSE)
            cnumber0 <- dim(sce)[2]
            sce <- Reduce(function(data, p_col) {
                idx <- substr(colData(data)[[p_col]], start = 1, stop = 1) == "T"
                data[, idx]
            }, produ_col, init = sce)
            cnumber1 <- dim(sce)[2]
            filtered <- cnumber0 - cnumber1
            message(sprintf("%d of cells filtered", filtered))
        }
    }
    ## retain only cells with allowed chain status
    if (!is.null(allowed_chain_status)) {
        message("checking allowed chain status...", appendLF = FALSE)
        cnumber0 <- dim(sce)[2]
        idx <- colData(sce)[["chain_status"]] %in% allowed_chain_status
        if (!any(idx)) {
            allowed_cs <- paste(allowed_chain_status, collapse = ", ")
            current_cs <- paste(unique(colData(sce)[["chain_status"]]), collapse = ", ")
            abort(sprintf(
                "Unsuitable allowed_chain_status,\n The current allowed_chain_status: %s.\n While the chain status in the dataset: %s.",
                allowed_cs, current_cs
            ))
        }
        sce <- sce[, idx]
        cnumber1 <- dim(sce)[2]
        filtered <- cnumber0 - cnumber1
        message(sprintf("%d of cells filtered", filtered))
    }
    ## subset sce by subsetby and groups
    if (!is.null(groups) && !is.null(subsetby)) {
        msg1 <- paste(as.character(substitute(groups))[-1], collapse = ", ")
        msg2 <- as.character(substitute(subsetby))
        message(sprintf("Subsetting data with %s in %s ...", msg1, msg2), appendLF = FALSE)
        cnumber0 <- dim(sce)[2]
        idx <- Reduce(`|`, lapply(groups, function(i) {
            colData(sce)[[subsetby]] %in%
                i
        }))
        sce <- sce[, idx]
        cnumber1 <- dim(sce)[2]
        filtered <- cnumber0 - cnumber1
        message(sprintf("%d of cells filtered", filtered))
    }
    ## extract main VDJ from specified columns
    message("VDJ data extraction begin:")
    if (is.null(extract_cols)) {
        message("Parameter extract_cols do not provided, automatically geneterate colnames for extraction.")
        if (!length(grep("_VDJ_main|_VJ_main", names(colData(sce))))) {
            v_call <- if ("v_call_genotyped_VDJ" %in% colnames(colData(sce))) {
                "v_call_genotyped_"
            } else {
                "v_call_"
            }
            prefix <- c(v_call, "d_call_", "j_call_")
            if (!is.null(mode_option)) {
                # can be pack as another function
                suffix <- c("_VDJ", "_VJ")
                extr_cols <- as.vector(outer(prefix, suffix, function(x, y) {
                    paste0(
                        x,
                        mode_option, y
                    )
                }))
                extr_cols <- extr_cols[extr_cols != paste0(
                    "d_call_", mode_option,
                    "_VJ"
                )]
            } else {
                suffix <- c("VDJ", "VJ")
                extr_cols <- as.vector(outer(prefix, suffix, function(x, y) {
                    paste0(
                        x,
                        y
                    )
                }))
                extr_cols <- extr_cols[extr_cols != paste0("d_call_", "VJ")]
            }
            msg <- paste(extr_cols, collapse = ", ")
            message(sprintf("Detect whether colData %s already exist...", msg))
            if (!any(extr_cols %in% colnames(colData(sce)))) {
                message(sprintf("Creating %s colData based on column CTgene", msg))
                splited_TCR <- splitCTgene(sce)
                if (length(splited_TCR[[1]]) != length(extr_cols)) {
                    abort(paste(
                        "Keyerror: Automatically generated colnames's length is",
                        length(extr_cols),
                        ".It must have the same number with the vdj data columns, which is",
                        length(splited_TCR[[1]]),
                        "\nYou could use parameter extract_cols to specify the columns to match the length"
                    ))
                } else {
                    vdj <- lapply(seq(length(extr_cols)), function(X, sc) {
                        vapply(X = sc, "[", X, FUN.VALUE = character(1))
                    }, sc = splited_TCR)
                    names(vdj) <- extr_cols
                    colData(sce) <- cbind(colData(sce), vdj)
                }
            } else if (!all(extr_cols %in% colnames(colData(sce)))) {
                abort(paste(
                    "Keyerror: Automatically generated colnames",
                    paste0(extr_cols[!extr_cols %in% colnames(colData(sce))],
                        collapse = ", "
                    ), "must be contained in colData",
                    "\nYou could use parameter extract_cols to specify the columns to extract TCR"
                ))
            }
            message(sprintf("Extract main TCR from %s ...", msg), appendLF = FALSE)
            sce <- Reduce(function(data, ex_col) {
                tem <- colData(data)[[ex_col]]
                strtem <- strsplit(as.character(tem), "\\|")
                colData(data)[[paste(ex_col, "main", sep = "_")]] <- vapply(strtem,
                    `[`, 1,
                    FUN.VALUE = character(1)
                )
                data
            }, extr_cols, init = sce)
            message("Complete.")
        }
    } else {
        msg <- paste(extract_cols, collapse = ", ")
        if (!any(extract_cols %in% colnames(colData(sce)))) {
            message(sprintf("ColData does not exist, Creating %s colData based on column CTgene", msg))
            splited_TCR <- splitCTgene(sce)
            if (length(splited_TCR[[1]]) != length(extract_cols)) {
                abort(paste(
                    "Keyerror: Colnames",
                    paste0(extract_cols[!extract_cols %in% colnames(colData(sce))],
                        collapse = ", "
                    ), "must have the same length with the vdj data, which is of the length",
                    length(splited_TCR[[1]]),
                    "\nYou could modify parameter extract_cols to specify the columns to match the length"
                ))
            } else {
                vdj <- lapply(seq(length(extract_cols)), function(X, sc) {
                    vapply(X = sc, "[", X, FUN.VALUE = character(1))
                }, sc = splited_TCR)
                names(vdj) <- extract_cols
                colData(sce) <- cbind(colData(sce), vdj)
            }
        } else if (!all(extract_cols %in% colnames(colData(sce)))) {
            abort(paste(
                "Keyerror: Colnames",
                paste0(extract_cols[!extract_cols %in% colnames(colData(sce))],
                    collapse = ", "
                ), "must be contained in colData",
                "\nYou could modify parameter extract_cols to specify the columns to extract TCR"
            ))
        }
        message(sprintf("Extract main TCR from %s ...", msg), appendLF = FALSE)
        sce <- Reduce(function(data, ex_col) {
            tem <- colData(data)[[ex_col]]
            strtem <- strsplit(as.character(tem), "\\|")
            colData(data)[[paste(ex_col, "main", sep = "_")]] <- vapply(strtem, `[`,
                1,
                FUN.VALUE = character(1)
            )
            data
        }, extract_cols, init = sce)
        message("Complete.")
    }
    # remove unclear mapping
    if (filter_unmapped) {
        filter_pattern = ",|None|No_contig"
        extr_cols <- c()
        if (!is.null(mode_option)) {
            if (any(check_vdj_mapping)) {
                vdj_mapping <- c('v_call', 'd_call', 'j_call')
                extr_cols <- c(extr_cols, paste(vdj_mapping[check_vdj_mapping], mode_option, "VDJ_main",
                    sep = "_"
                ))
            }
            if (any(check_vj_mapping)) {
                vj_mapping <- c('v_call', 'j_call')
                extr_cols <- c(extr_cols, paste(vj_mapping[check_vj_mapping], mode_option, "VJ_main",
                    sep = "_"
                ))
            }
        } else {
            if (is.null(extract_cols)) {
                if (any(check_vdj_mapping)) {
                    vdj_mapping <- c('v_call', 'd_call', 'j_call')
                    extr_cols <- c(extr_cols, paste(vdj_mappin[check_vdj_mapping], "VDJ_main",
                        sep = "_"
                    ))
                }
                if (any(check_vj_mapping)) {
                    vj_mapping <- c('v_call', 'j_call')
                    extr_cols <- c(extr_cols, paste(vj_mapping[check_vj_mapping], "VJ_main", sep = "_"))
                }
            } else {
                if (!is.null(check_extract_cols_mapping)) {
                    extr_cols <- check_extract_cols_mapping
                }
            }
        }
        if (!is.null(extr_cols)) {
            msg <- paste(extr_cols, collapse = ", ")
            message(sprintf("Filtering cells from %s ...", msg), appendLF = FALSE)
            cnumber0 <- dim(sce)[2]
            sce <- Reduce(function(x, y) {
                .filterCells(
                    sce = x, col_n = y, filter_pattern = filter_pattern,
                    remove_missing = remove_missing
                )
            }, extr_cols, init = sce)
            cnumber1 <- dim(sce)[2]
            filtered <- cnumber0 - cnumber1
            message(sprintf("%d of cells filtered", filtered))
        }
    }
    message(sprintf("%d of cells remain.", dim(sce)[2]))
    return(sce)
}
