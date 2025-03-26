# Project : ElsaBLab_fRagmentomics

#' Sanity Check for Mutation Input Data
#' This function verifies the validity of the mutation input data by
#' checking the chromosome, position, reference, and alternative alleles.
#' It filters out invalid rows and returns only the cleaned dataset.
#'
#' @param mut_info A data frame with mutation information.
#'
#' @return A filtered data frame containing only valid mutations.
#'
#' @keywords internal
sanity_check_read_mut <- function(mut_info) {
  # Initialize the cleaned dataframe
  read_mut_cleaned <- data.frame()

  for (i in seq_len(nrow(mut_info))) {
    chr <- mut_info[i, 1]
    pos <- mut_info[i, 2]
    ref <- mut_info[i, 3]
    alt <- mut_info[i, 4]

    # Validate inputs
    if (!check_chr_input(chr) || !check_pos_input(pos) || !check_ref_alt_input(ref, alt)) {
      warning(paste("Invalid row:", i, "- CHROM:", chr, "POS:", pos, "REF:", ref, "ALT:", alt))
      next
    }

    # Store valid rows
    read_mut_cleaned <- rbind(read_mut_cleaned, mut_info[i, , drop = FALSE])
  }

  # Check if final dataframe contains at least one row
  if (nrow(read_mut_cleaned) == 0) {
    stop("No valid mutations found after sanity check.")
  }

  read_mut_cleaned
}

#' Check Chromosome Input
#' It must be either "chrN" or "N" (where N is between 1 and 22).
#'
#' @param chr Character vector representing the chromosome of interest.
#'
#' @return TRUE if the chromosome is valid, otherwise FALSE.
#'
#' @noRd
check_chr_input <- function(chr) {
  if (is.null(chr) || is.na(chr) || chr == "") {
    return(FALSE)
  }

  if (grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y)$", chr) ||
    grepl("^([1-9]|1[0-9]|2[0-2]|X|Y)$", chr)) {
    return(TRUE)
  }

  return(FALSE)
}

#' Check Position Input
#' Ensures that the genomic position is a valid positive integer.
#'
#' @param pos Numeric value representing the Genomic position of interest.
#'
#' @return TRUE if the position is a valid integer greater than zero.
#'
#' @noRd
check_pos_input <- function(pos) {
  if (is.null(pos) || is.na(pos) || pos == "") {
    return(FALSE)
  }

  if (is.numeric(pos) && pos == as.integer(pos) && pos > 0) {
    return(TRUE)
  }

  return(FALSE)
}

#' Check Reference and Alternative Alleles
#' Ensures that the reference and alternative alleles follow specific rules
#'
#' @param ref Character vector representing reference base(s).
#' @param alt Character vector representing alternative base(s).
#'
#' @return TRUE if both alleles are valid, otherwise FALSE.
#'
#' @noRd
check_ref_alt_input <- function(ref, alt) {
  specific_values <- c("", ".", "-", "_", NA, "NA")

  if ((ref %in% specific_values) && (alt %in% specific_values)) {
    return(FALSE)
  }

  if (grepl(",", ref) || grepl(",", alt)) {
    return(FALSE)
  }

  valid_pattern <- "^[ATCG]*[._-]?$"
  if ((ref %in% specific_values || grepl(valid_pattern, ref)) &&
    (alt %in% specific_values || grepl(valid_pattern, alt))) {
    return(TRUE)
  }

  return(FALSE)
}
