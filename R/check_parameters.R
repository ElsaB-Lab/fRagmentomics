#' Check and validate all input parameters for the analysis pipeline
#'
#' @description This function serves as a centralized gatekeeper, validating all user-provided parameters before the main
#' analysis begins. It calls a series of specialized  helper functions to check each parameter for correctness (e.g., file existence,
#' data type, valid values) and stops execution with an informative error message if any check fails.
#'
#' @inheritParams run_fRagmentomics
#'
#' @return None. The function stops execution if error.
#'
#' @keywords internal
check_parameters <- function(
    mut,
    bam,
    fasta,
    sample_id,
    neg_offset_mate_search,
    pos_offset_mate_search,
    one_based,
    flag_bam_list,
    report_tlen,
    report_softclip,
    report_5p_3p_bases_fragment,
    remove_softclip,
    retain_fail_qc,
    tmp_folder,
    output_path,
    verbose,
    n_cores) {
  check_mut(mut)
  check_bam(bam, verbose)
  check_fasta(fasta, verbose)
  check_sample(sample_id)
  check_neg_offset_mate_search(neg_offset_mate_search)
  check_pos_offset_mate_search(pos_offset_mate_search)
  check_one_based(one_based)
  check_flag_bam_list(flag_bam_list)
  check_report_tlen(report_tlen)
  check_report_softclip(report_softclip)
  check_report_bases_fragm_5p_3p(report_5p_3p_bases_fragment)
  check_remove_softclip(remove_softclip)
  check_retain_fail_qc(retain_fail_qc)
  check_tmp_folder(tmp_folder)
  check_output_path(output_path)
  check_verbose(verbose)
  check_n_cores(n_cores)
}

#' Check if the mut file exists
#'
#' @inheritParams check_parameters
#'
#' @noRd
check_mut <- function(mut) {
  # Check that mut is a single value, whether a filepath or a str representation
  if (length(mut) != 1) {
    stop("The parameter 'mut' should be a single value, not multiple elements.")
  }

  # Check if mut is a valid file format (VCF, TSV, or their compressed versions)
  is_file_format <- grepl("\\.(vcf|tsv)(\\.gz)?$", mut)

  # If it's a file, check if it exists
  if (is_file_format && !file.exists(mut)) {
    stop("The Mutation file does not exist: ", mut)
  }
}

#' Check if the BAM file exists
#'
#' @inheritParams check_parameters
#'
#' @importFrom Rsamtools indexBam
#'
#' @noRd
check_bam <- function(bam, verbose) {
  # Check if the bam have a correct extension
  if (!grepl("\\.bam$", bam)) {
    stop("The file does not have a valid BAM extension (.bam): ", bam)
  }

  # Check if the BAM file exists
  if (!file.exists(bam)) {
    stop("The BAM file does not exist: ", bam)
  }

  # Define the expected BAM index file (.bai)
  bam_index <- paste0(bam, ".bai")

  # If the BAM index is missing, create it
  if (!file.exists(bam_index)) {
    if (verbose) {
      message("Creating BAM index...")
    }
    Rsamtools::indexBam(bam)
  }
}

#' Check if the FASTA file exists and has an index
#'
#' @inheritParams check_parameters
#'
#' @importFrom Rsamtools indexFa
#'
#' @noRd
check_fasta <- function(fasta, verbose) {
  if (!grepl("\\.fa(sta)?$", fasta)) {
    stop("The file does not have a valid FASTA extension (.fa or .fasta): ", fasta)
  }

  if (!file.exists(fasta)) {
    stop("The FASTA file does not exist: ", fasta)
  }

  fasta_index <- paste0(fasta, ".fai")
  if (!file.exists(fasta_index)) {
    if (verbose) {
      message("Creating FASTA index...")
    }
    Rsamtools::indexFa(fasta)
  }
}

#' Check if the sample ID is valid
#'
#' @inheritParams check_parameters
#'
#' @noRd
check_sample <- function(sample_id) {
  if (!is.na(sample_id) && sample_id == "") {
    stop("sample ID cannot be empty. Can be NA.")
  }
}

#' Check if the neg_offset_mate_search parameter is valid
#'
#' @inheritParams check_parameters
#'
#' @noRd
check_neg_offset_mate_search <- function(neg_offset_mate_search) {
  if (!is.integer(neg_offset_mate_search)) {
    stop("neg_offset_mate_search must be integer.")
  }
  if (neg_offset_mate_search > 0) {
    stop("neg_offset_mate_search must be negative or equal to 0.")
  }
}

#' Check if the pos_offset_mate_search parameter is valid
#'
#' @inheritParams check_parameters
#'
#' @noRd
check_pos_offset_mate_search <- function(pos_offset_mate_search) {
  if (!is.integer(pos_offset_mate_search)) {
    stop("pos_offset_mate_search must be interger.")
  }
  if (pos_offset_mate_search < 0) {
    stop("pos_offset_mate_search must be positive or equal to 0..")
  }
}

#' Check if the one_based parameter is valid
#'
#' @inheritParams check_parameters
#'
#' @noRd
check_one_based <- function(one_based) {
  if (!is.logical(one_based) || length(one_based) != 1) {
    stop("one_based must be a single logical value.")
  }
}

#' Check if the flag_bam_list parameter is valid
#'
#' @inheritParams check_parameters
#'
#' @importFrom Rsamtools scanBamFlag
#'
#' @noRd
check_flag_bam_list <- function(flag_bam_list) {
  # Check if the input is a list
  if (!is.list(flag_bam_list)) {
    stop("flag_bam_list' must be a list.")
  }

  # Check if all values in the list are logical (TRUE, FALSE, or NA)
  all_values_logical <- all(vapply(flag_bam_list, is.logical, logical(1)))
  if (!all_values_logical) {
    stop(" All values in 'flag_bam_list' must be logical (TRUE, FALSE, or NA).")
  }

  # Check if the list items are named and if the names are valid
  if (length(flag_bam_list) > 0) {
    valid_flag_names <- names(formals(Rsamtools::scanBamFlag))
    user_flag_names <- names(flag_bam_list)

    if (is.null(user_flag_names) || any(user_flag_names == "")) {
      stop("All elements in a non-empty 'flag_bam_list' must be named.")
    }

    invalid_names <- user_flag_names[!user_flag_names %in% valid_flag_names]
    if (length(invalid_names) > 0) {
      stop(sprintf(
        "Invalid name(s) found in 'flag_bam_list': %s.\n\nSee ?Rsamtools::scanBamFlag for a list of valid flag names.",
        paste(shQuote(invalid_names), collapse = ", ")
      ))
    }
  }
}

#' Check if the report_5p_bases parameter is valid
#'
#' @inheritParams check_parameters
#'
#' @noRd
check_report_bases_fragm_5p_3p <- function(report_5p_3p_bases_fragment) {
  if (!is.integer(report_5p_3p_bases_fragment)) {
    stop("report_bases_fragment_5p_3p must be integer.")
  }
  if (report_5p_3p_bases_fragment < 0) {
    stop("report_bases_fragment_5p_3p must be non-negative.")
  }
}

#' Check if the report_tlen parameter is valid
#'
#' @inheritParams check_parameters
#'
#' @noRd
check_report_tlen <- function(report_tlen) {
  if (!is.logical(report_tlen) || length(report_tlen) != 1) {
    stop("report_tlen must be a single logical value.")
  }
}

#' Check if the report_softclip parameter is valid
#'
#' @inheritParams check_parameters
#'
#' @noRd
check_report_softclip <- function(report_softclip) {
  if (!is.logical(report_softclip) || length(report_softclip) != 1) {
    stop("report_softclip must be a single logical value.")
  }
}

#' Check if the check_remove_softclip parameter is valid
#'
#' @inheritParams check_parameters
#'
#' @noRd
check_remove_softclip <- function(remove_softclip) {
  if (!is.logical(remove_softclip) || length(remove_softclip) != 1) {
    stop("remove_softclip must be a single logical value.")
  }
}


#' Check if the retain_fail_qc parameter is valid
#'
#' @inheritParams check_parameters
#'
#' @noRd
check_retain_fail_qc <- function(retain_fail_qc) {
  if (!is.logical(retain_fail_qc) || length(retain_fail_qc) != 1) {
    stop("retain_fail_qc must be a single logical value.")
  }
}

#' Check if the tmp folder exists. If not, create it.
#'
#' @inheritParams check_parameters
#'
#' @noRd
check_tmp_folder <- function(tmp_folder) {
  if (is.na(tmp_folder)) {
    return(tempdir())
  }

  if (!is.character(tmp_folder) || length(tmp_folder) != 1) {
    stop("tmp_folder must be a single character string or NA.")
  }

  if (!dir.exists(tmp_folder)) {
    dir.create(tmp_folder, recursive = TRUE, showWarnings = FALSE)
  }
}

#' Check if a directory exists. If not, create it.
#'
#' @param output_path A single character string for the directory path.
#'
#' @noRd
check_output_path <- function(output_path) {
  # Return silently if the path is empty, NA or NULL
  if (is.null(output_path) || is.na(output_path) || output_path == "") {
    return(invisible(NULL))
  }

  # Validate the input type
  if (!is.character(output_path) || length(output_path) != 1) {
    stop("output_path' must be a single character string.")
  }
}

#' Check if the verbose parameter is valid
#'
#' @inheritParams check_parameters
#'
#' @noRd
check_verbose <- function(verbose) {
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value.")
  }
}

#' Check if the n_cores parameter is valid
#'
#' @inheritParams check_parameters
#'
#' @noRd
check_n_cores <- function(n_cores) {
  if (!is.integer(n_cores)) {
    stop("n_cores must be integer.")
  }
  if (n_cores <= 0) {
    stop("n_cores must be positive.")
  }
}
