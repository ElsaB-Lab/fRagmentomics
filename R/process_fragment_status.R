# Project : ElsaBLab_fRagmentomics

#' Check if a read supports a mutation
#'
#' This helper function determines whether a given read contains a mutation, based on the type of mutation
#' and the read's base information. For SNVs, it checks for an exact match with the alternative allele.
#' For insertions and deletions, it checks whether the base string begins with "+" or "-", respectively.
#'
#' @param info A character string or NA representing the base or mutation detection status from a read.
#' @param mutation_type A character string indicating the type of mutation. Should be one of "SNV", "insertion", or "deletion".
#' @param alt A character string representing the alternative allele (used only for SNVs).
#'
#' @return A logical value: TRUE if the mutation is detected in the read, FALSE otherwise.
#'
#' @noRd
mutation_detected <- function(info, mutation_type, alt) {
  if (is.na(info)) {
    return(FALSE)
  }

  if (mutation_type == "SNV") {
    return(info == alt)
  } else if (mutation_type == "insertion") {
    # Check if the string starts with "+"
    return(grepl("^\\+", info))
  } else if (mutation_type == "deletion") {
    # Check if the string starts with "-"
    return(grepl("^\\-", info))
  } else {
    return(FALSE)
  }
}

#' Process Fragment Status
#'
#' This function determines the mutation status of a DNA fragment based on the mutation type and
#' read information. It supports Single Nucleotide Variants (SNV), insertions, and deletions.
#' The function returns a status based on whether both reads cover the mutation position or only one read covers it.
#'
#' @param alt A character string representing the alternative allele for SNV mutations.
#' @param mutation_type A character string indicating the type of mutation. Should be one of "SNV", "insertion", or "deletion".
#' @param r_info_base1 A character string or NA indicating the base or mutation detection status from the first read.
#' @param r_info_base2 A character string or NA indicating the base or mutation detection status from the second read.
#'
#' @return A character string representing the fragment status. Possible values are:
#'  -"MUT": Both reads support the mutation (or at least one read supports it when coverage is incomplete).
#'  -"WT_but_other_read_mut": Both reads are available and only one supports the mutation.
#'  -"WT": Neither read supports the mutation.
#'
#' @noRd
process_fragment_status <- function(alt, mutation_type, r_info_base1, r_info_base2) {
  # Check mutation detection for each read
  detected1 <- mutation_detected(r_info_base1, mutation_type, alt)
  detected2 <- mutation_detected(r_info_base2, mutation_type, alt)

  # Check whether both reads cover the position (i.e., are not NA)
  both_cover <- !is.na(r_info_base1) && !is.na(r_info_base2)

  if (both_cover) {
    if (detected1 && detected2) {
      fragment_status <- "MUT"
    } else if (detected1 || detected2) {
      fragment_status <- "WT_but_other_read_mut"
    } else {
      fragment_status <- "WT"
    }
  } else {
    if (detected1 || detected2) {
      fragment_status <- "MUT"
    } else {
      fragment_status <- "WT"
    }
  }

  return(fragment_status)
}
