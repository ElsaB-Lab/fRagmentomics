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
process_fragment_status <- function(alt, mutation_type, r_info_base1, r_info_base2){
  
  # Define target value based on mutation type
  target_value <- switch(mutation_type,
                         SNV = alt,
                         insertion = "insertion_detected",
                         deletion = "deletion_detected",
                         NULL)
                         
  if (is.null(target_value)) {
    return("Warning: Unknown mutation type")
  }
  
  # Check if both reads are present (i.e., not NA)
  both_cover <- !is.na(r_info_base1) && !is.na(r_info_base2)
  
  if (both_cover) {
    # When both reads are available:
    # If both reads match the target value -> MUT
    # If one matches -> WT_but_other_read_mut
    # Else -> WT
    if (r_info_base1 == target_value && r_info_base2 == target_value){
      fragment_status <- "MUT"
    } else if (r_info_base1 == target_value || r_info_base2 == target_value) {
      fragment_status <- "WT_but_other_read_mut"
    } else {
      fragment_status <- "WT"
    }
  } else {
    # When one or both reads are missing:
    # MUT if at least one read matches the target value; otherwise, WT
    if ((!is.na(r_info_base1) && r_info_base1 == target_value) || (!is.na(r_info_base2) && r_info_base2 == target_value)) {
      fragment_status <- "MUT"
    } else {
      fragment_status <- "WT"
    }
  }
  
  return(fragment_status)
}