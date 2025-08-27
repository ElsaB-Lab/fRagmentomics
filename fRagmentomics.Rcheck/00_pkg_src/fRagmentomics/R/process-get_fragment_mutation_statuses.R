#' Get detailed and simplified fragment mutation statuses.
#'
#' @description This function determines both the detailed and simplified mutation statuses of a DNA
#' fragment based on the mutation status of each of the two reads from the fragment.
#'
#' @param mstat_5p The mutation status of read 5p (e.g., "MUT:C>T_high_conf", "WT", "AMB_low_coverage").
#'                Can be "NA" for no coverage.
#' @param mstat_3p The mutation status of read 3p (e.g., "MUT:C>T_high_conf", "WT", "AMB_low_coverage").
#'                Can be "NA" for no coverage.
#'
#' @return A list containing two character strings:
#'         \itemize{
#'           \item 'Fragment_Status_Detail': A detailed status reflecting the combination,
#'                 retaining original text where applicable (e.g., "MUT & AMB_low_coverage").
#'           \item 'Fragment_Status_Simple': A simplified status (e.g., "MUT", "NON-TARGET MUT",
#'                 "DISCORDANT", "AMB", "ERROR").
#'         }
#'
#' @importFrom stats na.omit
#'
#' @keywords internal
get_fragment_mutation_statuses <- function(mstat_5p, mstat_3p) {
  # Internal helper to clean status for logical comparisons (removes extra detail)
  clean_status <- function(s) {
    if (is.na(s)) {
      return(NA_character_)
    }
    s_cleaned <- sub("(:.*|\\s.*)", "", s)
    return(s_cleaned)
  }

  base_mstat_5p <- clean_status(mstat_5p)
  base_mstat_3p <- clean_status(mstat_3p)

  # Initialize outputs
  fragment_status_detail <- NA_character_
  fragment_status_simple <- NA_character_

  # Define flags for cleaner logic
  is_na1 <- is.na(mstat_5p)
  is_na2 <- is.na(mstat_3p)

  is_mut1 <- base_mstat_5p == "MUT"
  is_mut2 <- base_mstat_3p == "MUT"
  is_wt1 <- base_mstat_5p == "WT"
  is_wt2 <- base_mstat_3p == "WT"
  is_amb1 <- base_mstat_5p == "AMB"
  is_amb2 <- base_mstat_3p == "AMB"
  is_other_mut1 <- base_mstat_5p == "OTH"
  is_other_mut2 <- base_mstat_3p == "OTH"

  # Function to combine original statuses for Fragment_Status_Detail
  # Handles sorting and unique values while preserving original string
  combine_original_statuses <- function(s1, s2) {
    original_statuses <- na.omit(c(s1, s2))
    if (length(original_statuses) == 0) {
      return(NA_character_)
    }
    if (length(original_statuses) == 1) {
      return(original_statuses[1])
    }
    if (original_statuses[1] == original_statuses[2]) {
      return(original_statuses[1])
    }

    # Sort based on the cleaned status, then get the original string
    sorted_pairs <- data.frame(
      clean = c(clean_status(s1), clean_status(s2)),
      original = c(s1, s2),
      stringsAsFactors = FALSE
    )
    sorted_pairs <- sorted_pairs[order(sorted_pairs$clean, sorted_pairs$original), ]

    # Filter out NA original strings and collapse unique ones
    unique_sorted_originals <- unique(na.omit(sorted_pairs$original))
    return(paste(unique_sorted_originals, collapse = " & "))
  }


  # 1. NA / NA
  if (is_na1 && is_na2) {
    fragment_status_detail <- "ERROR"
    fragment_status_simple <- "ERROR"
  }
  # 2. NA / MUT (and symmetrical)
  else if ((is_na1 && is_mut2) || (is_mut1 && is_na2)) {
    fragment_status_detail <- ifelse(is_na1, mstat_3p, mstat_5p)
    fragment_status_simple <- "MUT"
  }
  # 3. NA / WT (and symmetrical)
  else if ((is_na1 && is_wt2) || (is_wt1 && is_na2)) {
    fragment_status_detail <- ifelse(is_na1, mstat_3p, mstat_5p)
    fragment_status_simple <- "NON-TARGET MUT"
  }
  # 4. NA / AMB (and symmetrical)
  else if ((is_na1 && is_amb2) || (is_amb1 && is_na2)) {
    fragment_status_detail <- ifelse(is_na1, mstat_3p, mstat_5p)
    fragment_status_simple <- "AMB"
  }
  # 5. NA / OTH (and symmetrical)
  else if ((is_na1 && is_other_mut2) || (is_other_mut1 && is_na2)) {
    fragment_status_detail <- ifelse(is_na1, mstat_3p, mstat_5p)
    fragment_status_simple <- "NON-TARGET MUT"
  }
  # 6. MUT / MUT
  else if (is_mut1 && is_mut2) {
    fragment_status_detail <- combine_original_statuses(mstat_5p, mstat_3p)
    fragment_status_simple <- "MUT"
  }
  # 7. MUT / WT (and symmetrical) - DISCORDANT
  else if ((is_mut1 && is_wt2) || (is_wt1 && is_mut2)) {
    fragment_status_detail <- combine_original_statuses(mstat_5p, mstat_3p)
    fragment_status_simple <- "DISCORDANT"
  }
  # 8. MUT / AMB (and symmetrical) - Prioritize MUT
  else if ((is_mut1 && is_amb2) || (is_amb1 && is_mut2)) {
    fragment_status_detail <- combine_original_statuses(mstat_5p, mstat_3p)
    fragment_status_simple <- "MUT"
  }
  # 9. MUT / OTH (and symmetrical) - DISCORDANT
  else if ((is_mut1 && is_other_mut2) || (is_other_mut1 && is_mut2)) {
    fragment_status_detail <- combine_original_statuses(mstat_5p, mstat_3p)
    fragment_status_simple <- "DISCORDANT"
  }
  # 10. WT / WT
  else if (is_wt1 && is_wt2) {
    fragment_status_detail <- combine_original_statuses(mstat_5p, mstat_3p)
    fragment_status_simple <- "NON-TARGET MUT"
  }
  # 11. WT / AMB (and symmetrical) - Prioritize WT
  else if ((is_wt1 && is_amb2) || (is_amb1 && is_wt2)) {
    fragment_status_detail <- combine_original_statuses(mstat_5p, mstat_3p)
    fragment_status_simple <- "NON-TARGET MUT"
  }
  # 12. WT / OTH (and symmetrical) - DISCORDANT for detail, NON-TARGET MUT for simple
  else if ((is_wt1 && is_other_mut2) || (is_other_mut1 && is_wt2)) {
    fragment_status_detail <- combine_original_statuses(mstat_5p, mstat_3p)
    fragment_status_simple <- "NON-TARGET MUT"
  }
  # 13. AMB / AMB
  else if (is_amb1 && is_amb2) {
    fragment_status_detail <- combine_original_statuses(mstat_5p, mstat_3p)
    fragment_status_simple <- "AMB"
  }
  # 14. AMB / OTH (and symmetrical) - Prioritize OTH
  else if ((is_amb1 && is_other_mut2) || (is_other_mut1 && is_amb2)) {
    fragment_status_detail <- combine_original_statuses(mstat_5p, mstat_3p)
    fragment_status_simple <- "NON-TARGET MUT"
  }
  # 15. OTH / OTH
  else if (is_other_mut1 && is_other_mut2) {
    fragment_status_detail <- combine_original_statuses(mstat_5p, mstat_3p)
    fragment_status_simple <- "NON-TARGET MUT"
  }
  # Fallback for any unhandled or unexpected combination (should not happen)
  else {
    fragment_status_detail <- combine_original_statuses(mstat_5p, mstat_3p)
    fragment_status_simple <- "UNKNOWN"
  }

  # Return the results as a list
  return(list(
    Detail = fragment_status_detail,
    Simple = fragment_status_simple
  ))
}
