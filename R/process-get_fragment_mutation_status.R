#' Get fragment mutation status according to the mutation status of each read
#'
#' This function determines the mutation status of a DNA fragment based on
#' the mutation status of each read from the fragment.
#'
#' @param mstat_1 The mutation status of read 1
#' @param mstat_2 The mutation status of read 2
#'
#' @return A character string representing the fragment status:
#'
#' @keywords internal
get_fragment_mutation_status <- function(mstat_1, mstat_2){
  mstats <- c(mstat_1, mstat_2)

  if (all(mstats %in% c(NA))){
    "ERROR: position not covered"
  } else if (all(setequal(mstats, c("AMB", "WT")))){
    "WT"
  } else if (all(setequal(mstats, c("AMB", "MUT")))){
    "MUT"
  } else {
    paste(rev(sort(na.omit(mstats))), collapse = " & ")
  }
}
