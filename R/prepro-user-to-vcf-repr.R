# Project : ElsaBLab_fRagmentomics

#' Normalize Variant Data
#'
#' This function normalizes the REF and ALT columns by replacing invalid values with an empty string.
#' It also ensures that the CHROM column is treated as a string.
#'
#' @param df A data frame containing variant data with columns CHROM, POS, REF, and ALT.
#' @return A normalized data frame.
#' 
#' @noRd
normalize_variants <- function(df) {
  df$REF[df$REF %in% c("-", ".", "_", "NA")] <- ""
  df$ALT[df$ALT %in% c("-", ".", "_", "NA")] <- ""
  df$CHROM <- as.character(df$CHROM)
  return(df)
}

#' Normalize user-provided representation into VCF representation 
#'
#' This function normalizes the REF and ALT columns into the VCF format.
#' It also checks if all nucleotides in the REF column match those from a given FASTA reference.
#'
#' @param chr A string representing the chromosome (e.g., "1" or "chr1").
#' @param pos An integer representing the genomic position of the variant.
#' @param ref A string representing the reference allele at the given position.
#' @param alt A string representing the alternative allele without multiallelic sites.
#' @param fasta A reference genome in FASTA format, used to validate the REF nucleotide.
#'
#' @return A normalized VCF-like representation that can be used with BCT tools for further normalization.
#'
#' @noRd
normalize_user_rep_to_vcf_rep <- function(chr, pos, ref, alt, fasta) {
  
}
