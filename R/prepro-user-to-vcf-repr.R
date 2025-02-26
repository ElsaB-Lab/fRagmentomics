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