#' @title Retrieve mutation-specific information
#' @description Extracts base, quality, and indel information based on the mutation type.
#'
#' @param mutation_type Type of mutation: "mutation", "deletion", or "insertion".
#' @param pos Genomic position of interest.
#' @param ref Reference base.
#' @param alt Alternative base or sequence (for insertion).
#' @param del_info Deletion-specific information (if applicable).
#' @param read_stats A list of read-level statistics (position, cigar string, query sequence, quality scores).
#'
#' @return A list containing "base", "qual", and "indel" (if applicable).
#' 
#' @noRd 
get_mutation_info <- function(mutation_type, pos, ref, alt, del_info, read_stats) {
  
  if (mutation_type == "mutation") {
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
        del_info = del_info,
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
      base = NA,
      qual = NA,
      indel = "Error: mutation_type not defined properly"
    ))
  }
}
