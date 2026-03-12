#' Remove bad mutations
#'
#' @description This function verifies the validity of the mutation input data by checking the chromosome, position,
#' reference, and alternative alleles. It only returns rows that pass all checks.
#'
#' @param df_mut A dataframe with mutation information.
#'
#' @return A filtered dataframe containing only valid mutations.
#'
#' @keywords internal
remove_bad_mut <- function(df_mut) {
  # Build a logical vector of valid rows
  valid <- vapply(seq_len(nrow(df_mut)), function(i) {
    check_chr_input(df_mut[i, "CHROM"]) &&
      check_pos_input(df_mut[i, "POS"]) &&
      check_ref_alt_input(df_mut[i, "REF"], df_mut[i, "ALT"])
  }, logical(1))

  # Warn about invalid rows
  for (i in which(!valid)) {
    warning(sprintf(
      "Invalid row: %d - CHROM: %s POS: %s REF: %s ALT: %s",
      i, df_mut[i, "CHROM"], df_mut[i, "POS"], df_mut[i, "REF"], df_mut[i, "ALT"]
    ))
  }

  df_mut_clean <- df_mut[valid, , drop = FALSE]

  # Check if final dataframe contains at least one row
  if (nrow(df_mut_clean) == 0) {
    stop("No valid mutations found after sanity check.")
  }

  df_mut_clean
}

#' Check Chromosome Input
#' It must be either 'chrN' or 'N' (where N is between 1 and 22).
#'
#' @inheritParams process_fragment
#'
#' @return TRUE if the chromosome is valid, otherwise FALSE.
#'
#' @noRd
check_chr_input <- function(chr) {
  if (is.null(chr) || is.na(chr) || chr == "") {
    return(FALSE)
  }

  if (grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y)$", chr) || grepl(
    "^([1-9]|1[0-9]|2[0-2]|X|Y)$",
    chr
  )) {
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
#' @inheritParams process_fragment
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
  if ((ref %in% specific_values || grepl(valid_pattern, ref)) && (alt %in% specific_values ||
    grepl(valid_pattern, alt))) {
    return(TRUE)
  }

  return(FALSE)
}
