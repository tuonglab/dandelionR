#' setup_vdj_pseudobulk
#'
#' Function for data preprocessing. This function will filters the data based on productivity, chain status, subsets data, extracts main v(d)j, and removes the unmapping data
#' @param sce SingleCellExperiment object, vdj data should be contained in colData for filtering
#' @param mode_option optional, mode for extraction the V(D)J genes.
#' @param already.productive logical, whether the data is already filtered by productivity.
#' @param productive_cols character vector, names of colData used in productivity filtering, with NULL by default
#' @param productive_vj logical, with TRUE by default. Option in productivity filtering.  If True, cell will only be kept if the main VJ chain is productive
#' @param productive_vdj logical, with TRUE by default. Option in productivity filtering.  If True, cell will only be kept if the main VDJ chain is productive
#' @param allowed_chain_status character vectors, optional, if specified, the element should within c('single pair','Extra pair', 'Extra pair-exceptipn', 'Orphan VDJ','Orphan VDJ-exception')
#' @param subsetby character, with NULL by default, name of one colData provided for sub-setting.
#' @param groups character vector, with NULL by default, condition for sub-setting
#' @param extract_cols character vector, with NULL by default. colData names where VDJ/VJ information is stored, these colData will be used instead of the standard colData
#' @param filter_pattern character, optional, with ',|None|No_contig' by default.
#' @param check_vj_mapping character vector, optional. elements should come from c('v_call', 'j_call'), with c('v_call', 'j_call') by default.
#' @param check_vdj_mapping character vector, optional. elements should come from c('v_call','d_call', 'j_call'), with c('v_call', 'j_call') by default.
#' @param check_extract_cols_mapping character vecter, with NULL by default
#' @param remove_missing bool, True by default
#' @details
#' data will undergo several process, including productivity filtering, chain status filtering, subseting, main TCR extracting, and unmapping data removing.
#' - parameters for **productivity filtering**:
#'    - already.productive
#'      - If true, the function will skip the productivity filtering
#'      - If false, retain only productive entries in columns specified by mode_option or productive_cols
#'    - mode_option
#'      - function will check colData(s) named `productive_{mode_option}_{type}`, where type should be 'VDJ' or 'VJ' or both, depending on values of productive_vj and productive_vdj.
#'      - If set as 'NULL', the function needs the option 'extract_cols' to be specified
#'    - productive_cols
#'      - must be be specified when productivity filtering is need to conduct and mode_option is NULL.
#'      - where VDJ/VJ information is stored so that this will be used instead of the standard columns.
#'    - productive_vj, productive_vdj
#'      - If True, cell will only be kept if the main V(D)J chain is productive
#'  - parameter for **chain status filtering**: allowed_chain_status, chain status to be kept
#'  - parameters for **subsetting**: subsetby, groups
#'    - subsetting process will only be conducted when both parameters are provided. After subsetting, only the cell with \{groups\} feature in \{subsetby\} will be used for computing the VDJ feature space
#'  - parameter for **main v(d)j extraction**: extract_cols
#'  - parameters for **unmapping filtering**:
#'    - filter_pattern
#'      - pattern to be filtered form object.
#'      - If is set to be NULL, the unmaping filtering procees won't start
#'    - check_vj_mapping, check_vdj_mapping
#'      - only colData specified by these arguments (check_vj_mapping and check_vdj_mapping) will be checked for unclear mappings
#'    - check_extract_cols_mapping, related to extract_cols
#'      - Only colData specified by the argument will be checked for unclear mapping, the colData should first specified by extract_cols
#'    - remove_missing
#'      - If true, will remove cells with contigs matching the filter from the object.
#'      - If False, will mask them with a uniform value dependent on the column name.
#' @include check.R
#' @include filter.cells.R
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @return filtered SingleCellExperiment object
#' @export
setup_vdj_pseudobulk <- function(sce, mode_option = c("abT", "gdT", "B"), already.productive = TRUE,
    productive_cols = NULL, productive_vj = TRUE, productive_vdj = TRUE, allowed_chain_status = c("Single pair",
        "Extra pair", "Extra pair-exception", "Orphan VDJ", "Orphan VDJ-exception"),
    subsetby = NULL, groups = NULL, extract_cols = NULL, filter_pattern = ",|None|No_cotig",
    check_vj_mapping = c("v_call", "j_call"), check_vdj_mapping = c("v_call", "j_call"),
    check_extract_cols_mapping = NULL, remove_missing = TRUE) {
    # check if the data type is correct
    .class.check(sce, "SingleCellExperiment")
    mode_option <- match.arg(mode_option)
    .type.check(productive_cols, "character")
    .type.check(productive_vdj, "logical")
    .type.check(productive_vj, "logical")
    .type.check(subsetby, "character")
    .type.check(groups, "character")
    allowed_chain_status <- match.arg(allowed_chain_status, several.ok = TRUE)
    .type.check(extract_cols, "character")
    .type.check(filter_pattern, "character")
    check_vdj_mapping <- match.arg(check_vdj_mapping, c("v_call", "d_call", "j_call"),
        several.ok = TRUE)
    check_vj_mapping <- match.arg(check_vj_mapping, several.ok = TRUE)
    .type.check(check_extract_cols_mapping, "character")
    .type.check(remove_missing, "logical")

    # filtering retain only productive entries based on specified mode
    requireNamespace("rlang")
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
                rlang::abort("When mode_option is NULL, the productive_cols must be specified.")
            }
        } else {
            produ_col <- paste("productive", mode_option, c("VDJ", "VJ"), sep = "_")[c(productive_vdj,
                productive_vj)]
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
            rlang::abort(sprintf("Unsuitable allowed_chain_status,\n The current allowed_chain_status: %s.\n While the chain status in the dataset: %s.",
                allowed_cs, current_cs))
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
        idx <- Reduce(`|`, lapply(groups, function(i) colData(sce)[[subsetby]] %in%
            i))
        sce <- sce[, idx]
        cnumber1 <- dim(sce)[2]
        filtered <- cnumber0 - cnumber1
        message(sprintf("%d of cells filtered", filtered))
    }
    ## extract main VDJ from specified columns
    if (is.null(extract_cols)) {
        if (!length(grep("_VDJ_main|_VJ_main", names(colData(sce))))) {
            v_call <- if ("v_call_genotyped_VDJ" %in% colnames(colData(sce)))
                "v_call_genotyped_" else "v_call_"
            prefix <- c(v_call, "d_call_", "j_call_")
            if (!is.null(mode_option)) {
                # can be pack as another function
                suffix <- c("_VDJ", "_VJ")
                extr_cols <- as.vector(outer(prefix, suffix, function(x, y) paste0(x,
                  mode_option, y)))
                extr_cols <- extr_cols[extr_cols != paste0("d_call_", mode_option,
                  "_VJ")]
            } else {
                suffix <- c("VDJ", "VJ")
                extr_cols <- as.vector(outer(prefix, suffix, function(x, y) paste0(x,
                  y)))
                extr_cols <- extr_cols[extr_cols != paste0("d_call_", "VJ")]
            }
            msg <- paste(extr_cols, collapse = ", ")
            message(sprintf("Extract main TCR from %s ...", msg), appendLF = FALSE)
            sce <- Reduce(function(data, ex_col) {
                tem <- colData(data)[[ex_col]]
                strtem <- strsplit(as.character(tem), "\\|")
                colData(data)[[paste(ex_col, "main", sep = "_")]] <- vapply(strtem,
                  `[`, 1, FUN.VALUE = character(1))
                data
            }, extr_cols, init = sce)
            message("Complete.")
        }
    } else {
        msg <- paste(extract_cols, collapse = ", ")
        message(sprintf("Extract main TCR from %s ...", msg), appendLF = FALSE)
        sce <- Reduce(function(data, ex_col) {
            tem <- colData(data)[[ex_col]]
            strtem <- strsplit(as.character(tem), "\\|")
            colData(data)[[paste(ex_col, "main", sep = "_")]] <- vapply(strtem, `[`,
                1, FUN.VALUE = character(1))
            data
        }, extract_cols, init = sce)
        message("Complete.")
    }
    # remove unclear mapping
    if (!is.null(filter_pattern)) {
        extr_cols <- c()
        if (!is.null(mode_option)) {
            if (!is.null(check_vdj_mapping)) {
                extr_cols <- c(extr_cols, paste(check_vdj_mapping, mode_option, "VDJ_main",
                  sep = "_"))
            }
            if (!is.null(check_vj_mapping)) {
                extr_cols <- c(extr_cols, paste(check_vj_mapping, mode_option, "VJ_main",
                  sep = "_"))
            }
        } else {
            if (is.null(extract_cols)) {
                if (!is.null(check_vdj_mapping)) {
                  extr_cols <- c(extr_cols, paste(check_vdj_mapping, "VDJ_main",
                    sep = "_"))
                }
                if (!is.null(check_vj_mapping)) {
                  extr_cols <- c(extr_cols, paste(check_vj_mapping, "VJ_main", sep = "_"))
                }
            } else {
                if (!is.null(check_extract_cols_mapping))
                  extr_cols <- check_extract_cols_mapping
            }
        }
        if (!is.null(extr_cols)) {
            msg <- paste(extr_cols, collapse = ", ")
            message(sprintf("Filtering cells from %s ...", msg), appendLF = FALSE)
            cnumber0 <- dim(sce)[2]
            sce <- Reduce(function(x, y) {
                .filter.cells(sce = x, col_n = y, filter_pattern = filter_pattern,
                  remove_missing = remove_missing)
            }, extr_cols, init = sce)
            cnumber1 <- dim(sce)[2]
            filtered <- cnumber0 - cnumber1
            message(sprintf("%d of cells filtered", filtered))
        }
    }
    message(sprintf("%d of cells remain.", dim(sce)[2]))
    return(sce)
}

