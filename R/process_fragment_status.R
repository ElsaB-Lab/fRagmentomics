# Project : ElsaBLab_fRagmentomics

#' Check if a read supports a mutation
#' This helper function determines whether a given read contains a mutation,
#' based on the type of mutation and the read's base information.
#'
#' @inheritParams process_fragment
#' @param info A character string or NA indicating the base or
#'  mutation detection status
#'
#' @return TRUE if the mutation is detected in the read, FALSE otherwise.
#'
#' @noRd
mutation_detected <- function(ref, alt, info, mutation_type) {
  if (is.na(info)) {
    return(FALSE)
  }

  if (mutation_type == "SNV") {
    return(info == alt)
  } else if (mutation_type == "insertion") {
    # Check if info == alt for insertion
    return(substring(info, 2) == alt)
  } else if (mutation_type == "deletion") {
    # Check if info == ref for deletion
    return(substring(info, 2) == ref)
  } else {
    return(FALSE)
  }
}

#' Process Fragment Status
#' This function determines the mutation status of a DNA fragment based on
#' the mutation type and read information. It supports (SNV), insertions,
#' and deletions.
#'
#' @inheritParams process_fragment
#' @param r_info_base1 A character string or NA indicating the base or
#'  mutation detection status from the first read.
#' @param r_info_base2 A character string or NA indicating the base or
#'  mutation detection status from the second read.
#'
#' @return A character string representing the fragment status:
#'  -"MUT": Both reads support the mutation /
#'  (or at least one read supports it when coverage is incomplete).
#'  -"WT_but_other_read_mut": Both reads cover the position but only
#'  one supports the mutation.
#' -"WT_but_other_read_mut_with_other_alt": Both reads cover the position
#'  but one support a mutation that is not alt.
#' -"MUT_but_other_read_mut_with_other_alt": Both reads cover the position
#' but one is not ref or alt.
#' -"Error_both_read_mut_with_other_alt": Both reads cover the position
#' but both are not ref or alt.
#' -"Other_MUT": One read is covering and it's not ref or alt
#' -"WT": Neither read supports the mutation.
#'
#' @keywords internal
process_fragment_status <- function(ref,
                                    alt,
                                    mutation_type,
                                    r_info_base1,
                                    r_info_base2) {
  # If mutation_type is equal to deletion or insertion
  # Modify ref or alt to juste have the sequence of the indel
  if (mutation_type == "insertion") {
    alt <- substring(alt, 2)
  } else if (mutation_type == "deletion") {
    ref <- substring(ref, 2)
  }

  # Use the provided mutation_detected() function to check whether the read
  # supports the expected mutation
  mut_expected1 <- (!is.na(r_info_base1)) &&
    mutation_detected(ref, alt, r_info_base1, mutation_type)
  mut_expected2 <- (!is.na(r_info_base2)) &&
    mutation_detected(ref, alt, r_info_base2, mutation_type)

  # Helper to classify each read:
  # "mut"   : read shows the expected mutation (as defined by mutation_detected)
  # "wt"    : read equals the reference allele
  # "other" : read covers the position, but its value is neither ref nor alt
  get_status <- function(ref, mutation_type, r_info_base, mut_expected) {
    # Common check: if base is missing, return NA
    if (is.na(r_info_base)) {
      return(NA)
    }
    # If the mutation is expected, return "mut"
    if (mut_expected) {
      return("mut")
    }
    # For each mutation type, define when to return "wt"
    wt_condition <- switch(mutation_type,
      SNV = (r_info_base == ref),
      insertion = (r_info_base == "no_insertion_detected"),
      deletion = (r_info_base == "no_deletion_detected"),
      FALSE
    )
    if (wt_condition) {
      return("wt")
    }

    if (mutation_type == "insertion") {
      if (grepl("^\\+", r_info_base)) {
        return("Error. Insertion at the good position but different than alt")
      }
    } else if (mutation_type == "deletion") {
      if (grepl("^\\+", r_info_base)) {
        return("Error. Insertion at the good position but different than alt")
      }
    }

    # Fallback for any other case
    return("other")
  }

  status1 <- get_status(ref, mutation_type, r_info_base1, mut_expected1)
  status2 <- get_status(ref, mutation_type, r_info_base2, mut_expected2)

  # Determine if both reads cover the position
  both_cover <- !is.na(status1) && !is.na(status2)

  # If both reads cover the position:
  if (both_cover) {
    if (status1 == "mut" && status2 == "mut") {
      fragment_status <- "MUT"
    } else if ((status1 == "mut" && status2 == "wt") ||
      (status1 == "wt" && status2 == "mut")) {
      fragment_status <- "WT_but_other_read_mut"
    } else if ((status1 == "wt" && status2 == "other") ||
      (status1 == "other" && status2 == "wt")) {
      fragment_status <- "WT_but_other_read_mut_with_other_alt"
    } else if ((status1 == "mut" && status2 == "other") ||
      (status1 == "other" && status2 == "mut")) {
      fragment_status <- "MUT_but_other_read_mut_with_other_alt"
    } else if (status1 == "other" && status2 == "other") {
      fragment_status <- "Error_both_read_mut_with_other_alt"
    } else if (status1 == "wt" && status2 == "wt") {
      fragment_status <- "WT"
    } else {
      # Fallback for any unexpected combination
      fragment_status <- "Error_to_detect_ref/alt"
    }
  } else {
    # Incomplete coverage (only one read covers the position)
    if (is.na(status1) && is.na(status2)) {
      fragment_status <- "Error_1_read_should_cover_the_position"
    } else {
      # Incomplete coverage (only one read covers the position)
      if ((!is.na(status1) && status1 == "mut") ||
        (!is.na(status2) && status2 == "mut")) {
        fragment_status <- "MUT"
      } else if ((!is.na(status1) && status1 == "wt") ||
        (!is.na(status2) && status2 == "wt")) {
        fragment_status <- "WT"
      } else {
        fragment_status <- "Other_MUT"
      }
    }
  }

  return(fragment_status)
}
