#' Check input files for fRagmentomics function
#'
#' @inheritParams fRagmentomics
#'
#' @return None. The function stops execution if error.
#'
#' @keywords internal
check_input <- function(
    mut,
    bam,
    fasta,
    sample_id,
    neg_offset_mate_search,
    pos_offset_mate_search,
    one_based,
    flag_keep,
    flag_remove,
    report_tlen,
    report_softclip,
    report_5p_3p_bases_fragment,
    tmp_folder,
    n_cores) {
  check_mut(mut)
  check_bam(bam)
  check_fasta(fasta)
  check_sample(sample_id)
  check_neg_offset_mate_search(neg_offset_mate_search)
  check_pos_offset_mate_search(pos_offset_mate_search)
  check_one_based(one_based)
  check_flag_keep(flag_keep)
  check_flag_remove(flag_remove)
  check_report_tlen(report_tlen)
  check_report_softclip(report_softclip)
  check_report_bases_fragm_5p_3p(report_5p_3p_bases_fragment)
  check_tmp_folder(tmp_folder)
  check_n_cores(n_cores)
}

#' Check if the mut file exists
#'
#' @inheritParams check_input
#'
#' @noRd
check_mut <- function(mut) {
  # Check if mut is a valid file format (VCF, TSV, or their compressed versions)
  is_file_format <- grepl("\\.(vcf|tsv)(\\.gz)?$", mut)

  # If it's a file, check if it exists
  if (is_file_format && !file.exists(mut)) {
    stop("Error: The Mutation file does not exist: ", mut)
  }
}

#' Check if the BAM file exists
#'
#' @inheritParams check_input
#'
#' @importFrom Rsamtools indexBam
#'
#' @noRd
check_bam <- function(bam) {
  # Check if the BAM file exists
  if (!file.exists(bam)) {
    stop("Error: The BAM file does not exist: ", bam)
  }

  # Define the expected BAM index file (.bai)
  bam_index <- paste0(bam, ".bai")

  # If the BAM index is missing, create it
  if (!file.exists(bam_index)) {
    message("Creating BAM index...")
    Rsamtools::indexBam(bam)
  }
}

#' Check if the FASTA file exists and has an index
#'
#' @inheritParams check_input
#'
#' @importFrom Rsamtools indexFa
#'
#' @noRd
check_fasta <- function(fasta) {
  if (!file.exists(fasta)) {
    stop("Error: The FASTA file does not exist: ", fasta)
  }

  fasta_index <- paste0(fasta, ".fai")
  if (!file.exists(fasta_index)) {
    message("Creating FASTA index...")
    Rsamtools::indexFa(fasta)
  }
}

#' Check if the sample ID is valid
#'
#' @inheritParams check_input
#'
#' @noRd
check_sample <- function(sample_id) {
  if (!is.character(sample_id) && !is.na(sample_id)) {
    stop("Error: sample ID must be a character string or NA.")
  }
}

#' Check if the neg_offset_mate_search parameter is valid
#'
#' @inheritParams check_input
#'
#' @noRd
check_neg_offset_mate_search <- function(neg_offset_mate_search) {
  if (!is.integer(neg_offset_mate_search)) {
    stop("Error: neg_offset_mate_search must be integer.")
  }
  if (neg_offset_mate_search > 0) {
    stop("Error: neg_offset_mate_search must be negative or equal to 0.")
  }
}

#' Check if the pos_offset_mate_search parameter is valid
#'
#' @inheritParams check_input
#'
#' @noRd
check_pos_offset_mate_search <- function(pos_offset_mate_search) {
  if (!is.integer(pos_offset_mate_search)) {
    stop("Error: pos_offset_mate_search must be interger.")
  }
  if (pos_offset_mate_search < 0) {
    stop("Error: pos_offset_mate_search must be positive or equal to 0..")
  }
}

#' Check if the one_based parameter is valid
#'
#' @inheritParams check_input
#'
#' @noRd
check_one_based <- function(one_based) {
  if (!is.logical(one_based) || length(one_based) != 1) {
    stop("Error: one_based must be a single logical value.")
  }
}

#' Check if the flag_keep parameter is valid
#'
#' @inheritParams check_input
#'
#' @noRd
check_flag_keep <- function(flag_keep) {
  if (!is.numeric(flag_keep)) {
    stop("Error: flag_keep must be numeric.")
  }
  if (flag_keep < 0) {
    stop("Error: flag_keep must be non-negative.")
  }
}

#' Check if the flag_remove parameter is valid
#'
#' @inheritParams check_input
#'
#' @noRd
check_flag_remove <- function(flag_remove) {
  if (!is.numeric(flag_remove)) {
    stop("Error: flag_remove must be numeric.")
  }
  if (flag_remove < 0) {
    stop("Error: flag_remove must be non-negative.")
  }
}

#' Check if the report_5p_bases parameter is valid
#'
#' @inheritParams check_input
#'
#' @noRd
check_report_bases_fragm_5p_3p <- function(report_5p_3p_bases_fragment) {
  if (!is.integer(report_5p_3p_bases_fragment)) {
    stop("Error: report_bases_fragment_5p_3p must be integer.")
  }
  if (report_5p_3p_bases_fragment < 0) {
    stop("Error: report_bases_fragment_5p_3p must be non-negative.")
  }
}

#' Check if the report_tlen parameter is valid
#'
#' @inheritParams check_input
#'
#' @noRd
check_report_tlen <- function(report_tlen) {
  if (!is.logical(report_tlen) || length(report_tlen) != 1) {
    stop("Error: report_tlen must be a single logical value.")
  }
}

#' Check if the report_softclip parameter is valid
#'
#' @inheritParams check_input
#'
#' @noRd
check_report_softclip <- function(report_softclip) {
  if (!is.logical(report_softclip) || length(report_softclip) != 1) {
    stop("Error: report_softclip must be a single logical value.")
  }
}

#' Check if the tmp folder exists. If not, create it.
#'
#' @inheritParams check_input
#'
#' @noRd
check_tmp_folder <- function(tmp_folder) {
  if (is.na(tmp_folder)) {
    return(tempdir())
  }

  if (!is.character(tmp_folder) || length(tmp_folder) != 1) {
    stop("Error: tmp_folder must be a single character string or NA.")
  }

  if (!dir.exists(tmp_folder)) {
    dir.create(tmp_folder, recursive = TRUE, showWarnings = FALSE)
  }

  tmp_folder
}

#' Check if the n_cores parameter is valid
#'
#' @inheritParams check_input
#'
#' @noRd
check_n_cores <- function(n_cores) {
  if (!is.integer(n_cores)) {
    stop("Error: n_cores must be integer.")
  }
  if (n_cores <= 0) {
    stop("Error: n_cores must be positive.")
  }
}
