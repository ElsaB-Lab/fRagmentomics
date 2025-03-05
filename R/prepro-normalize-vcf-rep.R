# Project : ElsaBLab_fRagmentomics

#' Create a temporary VCF file
#'
#' This function creates a temporary VCF file from given variant information.
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
    temp_vcf <- tempfile(fileext = ".vcf")

    # Define VCF header
    vcf_header <- "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

    # Define VCF content
    vcf_content <- paste(chr, pos, ".", ref, alt, ".", ".", ".", sep = "\t")

    # Write to temporary file
    writeLines(c(vcf_header, vcf_content), temp_vcf)

    return(temp_vcf)
}


normalize_vcf_rep <- function(chr, pos, ref, alt, fasta) {
    
}