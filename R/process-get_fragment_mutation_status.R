#' Get fragment mutation status according to the mutation status of each read.
#'
#' This function determines the mutation status of a DNA fragment based on
#' the mutation status of each read from the fragment.
#'
#' @param fragment_status_broad The broad label mutation status of the fragment
#'
#' @return A character string representing the fragment status:
#'
#' @keywords internal
get_fragment_mutation_status <- function(fragment_status_broad) {
  if (fragment_status_broad == "ERROR: position not covered") {
    "ERROR: position not covered"
  } else if (fragment_status_broad == "WT" || fragment_status_broad == "Other MUT" || fragment_status_broad == "WT & Other MUT") {
    "Non-target MUT"
  } else if (fragment_status_broad == "MUT") {
    "MUT"
  } else if (fragment_status_broad == "WT & MUT" || fragment_status_broad == "Other MUT & MUT") {
    "DISCORDANT"
  } else {
    "AMB"
  }
}

#' Get fragment mutation status according to the mutation status of each read.
#' Detailed label of get_fragment_mutation_status
#'
#' This function determines the mutation status of a DNA fragment based on
#' the mutation status of each read from the fragment.
#'
#' @param mstat_1 The mutation status of read 1
#' @param mstat_2 The mutation status of read 2
#'
#' @return A character string representing the fragment status:
#'
#' @importFrom stats na.omit
#'
#' @keywords internal
get_fragment_mutation_status_broad <- function(mstat_1, mstat_2) {
  mstats <- c(mstat_1, mstat_2)

  if (all(mstats %in% c(NA))) {
    "ERROR: position not covered"
  } else if (all(setequal(mstats, c("AMB", "WT")))) {
    "WT"
  } else if (all(setequal(mstats, c("AMB", "MUT")))) {
    "MUT"
  } else if (all(setequal(mstats, c("AMB", "Other MUT")))) {
    "Other MUT"
  } else {
    paste(unique(rev(sort(na.omit(mstats)))), collapse = " & ")
  }
}
