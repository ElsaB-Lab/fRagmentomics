# Project : ElsaBLab_fRagmentomics

#' @title Retrieve mutation-specific information
#' @description Extracts base, quality, and indel information based on the mutation type.
#'
#' @param mutation_type Type of mutation: "SNV", "deletion", or "insertion".
#' @param pos Genomic position of interest.
#' @param ref Reference base.
#' @param alt Alternative base or sequence (for insertion).
#' @param read_stats A list of read-level statistics (position, cigar string, query sequence, quality scores).
#'
#' @return A list containing "base", "qual".
#' Base could be:
#'  - NA if the read doesn't cover the position of interest
#'  - Nucleodide for SNV if the read covers the position of interest
#'  - "insertion_detected" or "no_insertion_detected"
#'  - "deletion_detected" or "no_deletion_detected"
#' 
#' Qual could be:
#'  - NA if the read doesn't cover the position of interest
#'  - Qual for SNV and insertion if the read covers the position of interest
#'  - "no_insertion_detected"
#'  - "-" if the read covers the position of the deletion of interest or "no_deletion_detected"
#' 
#' @noRd 
get_mutation_info <- function(mutation_type, pos, ref, alt, read_stats) {
  
  if (mutation_type == "SNV") {
    return(
      get_base_qual_from_read(
        pos      = pos,
        r_pos    = read_stats$pos,
        r_cigar  = read_stats$cigar,
        r_query  = read_stats$query,
        r_qual   = read_stats$qual
      )
    )
    
  } else if (mutation_type == "deletion") {
    return(
      get_deletion(
        pos      = pos,
        ref      = ref,
        r_pos    = read_stats$pos,
        r_cigar  = read_stats$cigar
      )
    )
    
  } else if (mutation_type == "insertion") {
    return(
      get_insertion(
        pos      = pos,
        alt      = alt,
        r_pos    = read_stats$pos,
        r_cigar  = read_stats$cigar,
        r_query  = read_stats$query,
        r_qual   = read_stats$qual
      )
    )
    
  } else {
    return(list(
      base = "Error: mutation_type not defined properly",
      qual = "Error: mutation_type not defined properly"
    ))
  }
}
