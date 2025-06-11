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
#' the mutation type and read information. It supports single nucleotide variants (SNV),
#' insertions, and deletions, including an "ambiguous" status for indels.
#'
#' @inheritParams process_fragment
#' @param r_info_base1 A character string or NA indicating the base or
#'  mutation detection status from the first read. Can be "ambiguous".
#' @param r_info_base2 A character string or NA indicating the base or
#'  mutation detection status from the second read. Can be "ambiguous".
#'
#' @return A character string representing the fragment status:
#'  -"MUT": Both reads support the mutation, or one read is MUT and the other is Ambiguous or NA.
#'  -"WT": Both reads are wild-type, or one is WT and the other is Ambiguous or NA.
#'  -"Ambiguous": Both reads are ambiguous, or one is ambiguous and the other is NA.
#'  -"WT_but_other_read_mut": Both reads cover the position but only one supports the mutation.
#'  -"WT_but_other_read_mut_with_other_alt": Both reads cover the position but one supports a mutation that is not alt.
#'  -"MUT_but_other_read_mut_with_other_alt": Both reads cover the position but one is not ref or alt.
#'  -"Error_both_read_mut_with_other_alt": Both reads cover the position but both are not ref or alt.
#'  -"Other_MUT": One read is covering and it's not ref or alt.
#'
#' @keywords internal
process_fragment_status <- function(ref,
                                    alt,
                                    mutation_type,
                                    r_info_base1,
                                    r_info_base2) {
  # If mutation_type is equal to deletion or insertion
  # Modify ref or alt to just have the sequence of the indel
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
  # "mut"      : read shows the expected mutation (as defined by mutation_detected)
  # "wt"       : read equals the reference allele
  # "ambiguous": read has an ambiguous indel status
  # "other"    : read covers the position, but its value is neither ref nor alt
  get_status <- function(ref, mutation_type, r_info_base, mut_expected) {
    # Common check: if base is missing, return NA
    if (is.na(r_info_base)) {
      return(NA)
    }

    # Check for "ambiguous" status
    if (r_info_base == "ambiguous") {
      return("ambiguous")
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

    # Fallback for any other case
    return("other")
  }

  status1 <- get_status(ref, mutation_type, r_info_base1, mut_expected1)
  status2 <- get_status(ref, mutation_type, r_info_base2, mut_expected2)

  # Determine if both reads cover the position
  both_cover <- !is.na(status1) && !is.na(status2)

  # If both reads cover the position:
  if (both_cover) {
    # Sort statuses to handle symmetric cases (e.g., "wt" + "mut" is same as "mut" + "wt")
    statuses <- sort(c(status1, status2))
    s1 <- statuses[1]
    s2 <- statuses[2]

    # Handle "ambiguous" cases first
    if (s1 == "ambiguous") {
      if (s2 == "ambiguous") {
        fragment_status <- "Ambiguous" # Ambiguous + Ambiguous
      } else if (s2 == "mut") {
        fragment_status <- "MUT" # Ambiguous + MUT
      } else if (s2 == "wt") {
        fragment_status <- "WT" # Ambiguous + WT
      } else {
        # Fallback for ambiguous + other
        fragment_status <- "Error_to_detect_ref/alt"
      }
    } else if (s1 == "mut" && s2 == "mut") {
      fragment_status <- "MUT"
    } else if (s1 == "mut" && s2 == "wt") {
      fragment_status <- "WT_but_other_read_mut"
    } else if (s1 == "other" && s2 == "wt") {
      fragment_status <- "WT_but_other_read_mut_with_other_alt"
    } else if (s1 == "mut" && s2 == "other") {
      fragment_status <- "MUT_but_other_read_mut_with_other_alt"
    } else if (s1 == "other" && s2 == "other") {
      fragment_status <- "Error_both_read_mut_with_other_alt"
    } else if (s1 == "wt" && s2 == "wt") {
      fragment_status <- "WT"
    } else {
      # Fallback for any unexpected combination
      fragment_status <- "Error_to_detect_ref/alt"
    }
  } else {
    # Incomplete coverage (only one read covers the position)
    single_status <- if (!is.na(status1)) status1 else status2

    if (is.na(single_status)) {
      fragment_status <- "Error_1_read_should_cover_the_position"
    } else {
      # Add "ambiguous" case for single read coverage
      fragment_status <- switch(single_status,
        "mut" = "MUT",
        "wt" = "WT",
        "ambiguous" = "Ambiguous", # Ambiguous + NA
        "Other_MUT" # Default for "other"
      )
    }
  }

  return(fragment_status)
}
