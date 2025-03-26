#' Check input files for fRagmentomics function
#'
#' @inheritParams fRagmentomics
#' @param mut Could be a .vcf, .tsv or chr:pos:ref:alt.
#' @param bam A dataframe containing sequencing reads.
#' @param fasta Character vector.
#' @param sample_id Sample identifier.
#' @param neg_offset_mate_search Integer. Use in read_bam.
#'  Represents the number of nucleotides to extend upstream (negative direction)
#'  from the position of interest when querying the BAM file with Rsamtools.
#'  his extension ensures that paired reads are retrieved, even if only one mate
#'  overlaps the queried position.
#' @param pos_offset_mate_search Integer. Use in read_bam.
#' @param one_based Boolean. TRUE if fasta is in one based. False if in 0 based.
#' @param flag_keep Character vector. Use in read_bam.
#'  Represents the SAM flags that should be kept when filtering alignments.
#' @param flag_remove Character vector. Use in read_bam.
#'  Represents the SAM flags that should be excluded when filtering alignments.
#' @param report_tlen Boolean. Whether to include the TLEN (template length)
#'  information in the output.
#' @param report_softclip Boolean. Whether to include the number of soft-clipped
#'  bases at the fragment extremities in the output.
#' @param report_5p_3p_bases_fragment Integer. Whether to include N fragment
#'  extremity bases in the output.
#' @param n_cores Number of cores for parallel computation.
#'
#' @return None. The function stops execution if error.
#'
#' @keywords internal
check_input <- function(
    mut, bam, fasta, sample_id, neg_offset_mate_search,
    pos_offset_mate_search, one_based, flag_keep, flag_remove,
    report_tlen, report_softclip, report_5p_3p_bases_fragment, n_cores) {
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
  check_n_cores(n_cores)
}

#' Check if the mut file exists
#'
#' @inheritParams check_input
#' @param mut Could be a .vcf, .tsv or chr:pos:ref:alt.
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
#' @param bam Character. Path to the BAM file.
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
#' @param fasta Character. Path to the FASTA file.
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
#' @param sample_id Character or NA. Sample ID to be checked.
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
#' @param neg_offset_mate_search Numeric. The negative offset mate search value.
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
#' @param pos_offset_mate_search Numeric. The positive offset mate search value.
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
#' @param one_based Logical. Indicates whether indexing is one-based.
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
#' @param flag_keep Numeric. The flag to keep.
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
#' @param flag_remove Numeric. The flag to remove.
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
#' @param report_bases_fragment_5p_3p Integer. Whether to include N fragment
#'  extremity bases in the output.
#'
#' @noRd
check_report_bases_fragm_5p_3p <- function(report_5p_3p_bases_fragment) {
  if (!is.numeric(report_5p_3p_bases_fragment)) {
    stop("Error: report_bases_fragment_5p_3p must be numeric.")
  }
  if (report_5p_3p_bases_fragment < 0) {
    stop("Error: report_bases_fragment_5p_3p must be non-negative.")
  }
}

#' Check if the report_tlen parameter is valid
#'
#' @inheritParams check_input
#' @param report_tlen Logical. Indicates whether to report TLEN.
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
#' @param report_softclip Boolean. Whether to include the number of soft-clipped
#'  bases at the fragment extremities in the output.
#'
#' @noRd
check_report_softclip <- function(report_softclip) {
  if (!is.logical(report_softclip) || length(report_softclip) != 1) {
    stop("Error: report_softclip must be a single logical value.")
  }
}

#' Check if the n_cores parameter is valid
#'
#' @inheritParams check_input
#' @param n_cores Numeric. The number of cores to use.
#'
#' @noRd
check_n_cores <- function(n_cores) {
  if (!is.numeric(n_cores)) {
    stop("Error: n_cores must be numeric.")
  }
  if (n_cores <= 0) {
    stop("Error: n_cores must be positive.")
  }
}
