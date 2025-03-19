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
check_input <- function(mut, bam, fasta) {
  check_mut(mut)
  check_bam(bam)
  check_fasta(fasta)
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
