#' Check input files for fRagmentomics function
#'
#' This function verifies the existence of the BAM and FASTA files and ensures the FASTA file has an index.
#'
#' @param bam Character. Path to the BAM file.
#' @param fasta Character. Path to the FASTA file.
#'
#' @return None. The function stops execution if files are missing or creates an index for the FASTA file if needed.
#' 
#' @noRd 
check_input <- function(mut, bam, fasta, sample, neg_offset_mate_search, pos_offset_mate_search,
                        one_based, flag_keep, flag_remove, report_tlen, report_softclip, 
                        report_bases_fragment_5p_3p, n_cores) {
  check_mut(mut)
  check_bam(bam)
  check_fasta(fasta)
  check_sample(sample)
  check_neg_offset_mate_search(neg_offset_mate_search)
  check_pos_offset_mate_search(pos_offset_mate_search)
  check_one_based(one_based)
  check_flag_keep(flag_keep)
  check_flag_remove(flag_remove)
  check_report_tlen(report_tlen)
  check_report_softclip(report_softclip)
  check_report_bases_fragment_5p_3p(report_bases_fragment_5p_3p)
  check_report_bases_fragment_5p_3p(report_bases_fragment_5p_3p)
  check_n_cores(n_cores)
}

#' Check if the mut file exists
#'
#' This function verifies whether the specified mutation file exists. If not, it stops execution with an error message.
#'
#' @param bam Character. Path to the mutation file.
#'
#' @return None. The function stops execution if the file is missing.
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
#' This function verifies whether the specified BAM file exists. If not, it stops execution with an error message.
#'
#' @param bam Character. Path to the BAM file.
#'
#' @return None. The function stops execution if the file is missing.
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
#' This function verifies whether the specified FASTA file exists. If an index (.fai) is missing, it generates one using the Rsamtools package.
#'
#' @param fasta Character. Path to the FASTA file.
#'
#' @return None. The function stops execution if the file is missing or generates an index if needed.
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
#' This function verifies whether the sample ID is a character string or NA.
#' If not, it stops execution with an error message.
#'
#' @param sample Character or NA. Sample ID to be checked.
#'
#' @return None. The function stops execution if the sample ID is invalid.
#'
#' @noRd
check_sample <- function(sample) {
    if (!is.character(sample) && !is.na(sample)) {
        stop("Error: sample ID must be a character string or NA.")
    }
}

#' Check if the neg_offset_mate_search parameter is valid
#'
#' This function verifies that the provided negative offset mate search value is numeric and negative.
#'
#' @param neg_offset_mate_search Numeric. The negative offset mate search value.
#'
#' @return None. The function stops execution if the parameter is invalid.
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
#' This function verifies that the provided positive offset mate search value is numeric and positive.
#'
#' @param pos_offset_mate_search Numeric. The positive offset mate search value.
#'
#' @return None. The function stops execution if the parameter is invalid.
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
#' This function verifies that the provided one_based parameter is a single logical value.
#'
#' @param one_based Logical. Indicates whether indexing is one-based.
#'
#' @return None. The function stops execution if the parameter is invalid.
#'
#' @noRd
check_one_based <- function(one_based) {
  if (!is.logical(one_based) || length(one_based) != 1) {
    stop("Error: one_based must be a single logical value.")
  }
}

#' Check if the flag_keep parameter is valid
#'
#' This function verifies that the provided flag_keep parameter is numeric and non-negative.
#'
#' @param flag_keep Numeric. The flag to keep.
#'
#' @return None. The function stops execution if the parameter is invalid.
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
#' This function verifies that the provided flag_remove parameter is numeric and non-negative.
#'
#' @param flag_remove Numeric. The flag to remove.
#'
#' @return None. The function stops execution if the parameter is invalid.
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
#' This function verifies that the provided report_5p_bases parameter is numeric and non-negative.
#'
#' @param report_bases_fragment_5p_3p Numeric. The number of bases to report for the 5' read.
#'
#' @return None. The function stops execution if the parameter is invalid.
#'
#' @noRd
check_report_bases_fragment_5p_3p <- function(report_bases_fragment_5p_3p) {
  if (!is.numeric(report_bases_fragment_5p_3p)) {
    stop("Error: report_bases_fragment_5p_3p must be numeric.")
  }
  if (report_bases_fragment_5p_3p < 0) {
    stop("Error: report_bases_fragment_5p_3p must be non-negative.")
  }
}

#' Check if the report_tlen parameter is valid
#'
#' This function verifies that the provided report_tlen parameter is a single logical value.
#'
#' @param report_tlen Logical. Indicates whether to report TLEN.
#'
#' @return None. The function stops execution if the parameter is invalid.
#'
#' @noRd
check_report_tlen <- function(report_tlen) {
  if (!is.logical(report_tlen) || length(report_tlen) != 1) {
    stop("Error: report_tlen must be a single logical value.")
  }
}

#' Check if the report_softclip parameter is valid
#'
#' This function verifies that the provided report_softclip parameter is a single logical value.
#'
#' @param report_softclip Logical. Indicates whether to report softclip information.
#'
#' @return None. The function stops execution if the parameter is invalid.
#'
#' @noRd
check_report_softclip <- function(report_softclip) {
  if (!is.logical(report_softclip) || length(report_softclip) != 1) {
    stop("Error: report_softclip must be a single logical value.")
  }
}

#' Check if the n_cores parameter is valid
#'
#' This function verifies that the provided n_cores parameter is numeric and positive.
#'
#' @param n_cores Numeric. The number of cores to use.
#'
#' @return None. The function stops execution if the parameter is invalid.
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
