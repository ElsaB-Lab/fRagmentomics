# Project: ElsaBLab_fRagmentomics

#' Create a temporary VCF file
#'
#' This function creates a temporary VCF file from the given variant information.
#'
#' @param chr Character. The chromosome identifier.
#' @param pos Integer. The position of the variant.
#' @param ref Character. The reference allele.
#' @param alt Character. The alternative allele.
#'
#' @return Character. The path to the temporary VCF file.
#'
#' @noRd
create_temporary_vcf <- function(chr, pos, ref, alt) {
  # Create a temporary file with .vcf extension
  tmp_vcf <- tempfile(fileext = ".vcf")
  
  # Define VCF header
  vcf_header <- "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  
  # Define VCF content (one variant per line)
  vcf_content <- paste(chr, pos, ".", ref, alt, ".", ".", ".", sep = "\t")
  
  # Write the VCF header and content to the temporary file
  writeLines(c(vcf_header, vcf_content), tmp_vcf)
  
  return(tmp_vcf)
}

#' Normalize VCF file using bcftools norm
#'
#' This function normalizes a VCF file by calling bcftools norm on a temporary VCF file.
#' The normalized output is then read into a dataframe, extracting the columns:
#' chr, pos, ref, and alt. This function also handles complex cases including multi-allelic variants,
#' and it automatically deletes temporary files created during the process.
#'
#' @param chr Character. The chromosome identifier.
#' @param pos Integer. The position of the variant.
#' @param ref Character. The reference allele.
#' @param alt Character. The alternative allele.
#' @param fasta Character. The path to the reference FASTA file.
#'
#' @return A dataframe containing normalized variant information with columns: chr, pos, ref, alt.
#'
#' @noRd
normalize_vcf_rep <- function(chr, pos, ref, alt, fasta) {
    # Create a temporary VCF file from the input variant information
    tmp_vcf <- create_temporary_vcf(chr, pos, ref, alt)

    print("DEBUG: Temporary VCF File Content")
    print(readLines(tmp_vcf))

    # Create a temporary file to store the output of bcftools norm
    tmp_out_vcf <- tempfile(fileext = ".vcf")

    # Ensure that temporary files are removed after the function executes
    on.exit({
    if (file.exists(tmp_vcf)) file.remove(tmp_vcf)
    if (file.exists(tmp_out_vcf)) file.remove(tmp_out_vcf)
    }, add = TRUE)

    # Build the bcftools norm command. The -m +both option splits multi-allelic sites.
    cmd <- sprintf("bcftools norm -m +both -d all --check REF,ATL -f %s -o %s %s",
                shQuote(fasta),
                shQuote(tmp_out_vcf),
                shQuote(tmp_vcf))

    print(paste("Executing command:", cmd))
    system("bcftools --version", intern = TRUE)

    # Execute the command to normalize the VCF file
    system(cmd, intern = FALSE, ignore.stderr = TRUE)

    # Read the normalized VCF output into a dataframe.
    # It skips header lines (starting with "#") and selects relevant columns.
    normalized_variants <- readr::read_tsv(tmp_out_vcf, comment = "#", col_names = FALSE) %>%
        dplyr::select(X1, X2, X4, X5) %>%  # Selecting CHROM, POS, REF, ALT columns
        dplyr::rename(chr = X1, pos = X2, ref = X4, alt = X5)

    print("DEBUG: Normalized Variants Dataframe")
    print(normalized_variants)

    # Return the dataframe with normalized variant data
    return(normalized_variants)
}
