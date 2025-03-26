# Project : ElsaBLab_fRagmentomics

#' Retrieve mutation-specific information:
#' Extracts base, quality, and indel informations based on the mutation type.
#'
#' @inheritParams process_fragment
#' @param mutation_type In "SNV", "deletion", or "insertion".
#' @param pos Numeric value representing the Genomic position of interest.
#' @param ref Character vector representing reference base(s).
#' @param alt Character vector representing alternative base(s).
#' @param read_stats A list of read-level statistics.
#'
#' @return A list containing "base", "qual".
#' Base could be:
#'  - NA if the read doesn't cover the position of interest
#'  - Nucleodide for SNV if the read covers the position of interest
#'  - The sequence of the deletion or insertion if detected
#'  - "no_deletion_detected" or "no_insertion_detected"
#'
#' Qual could be:
#'  - NA if the read doesn't cover the position of interest
#'  - Qual for SNV if the read covers the position of interest
#'  - Qual of nucleotide before the indel if detected
#'  - "no_deletion_detected" or "no_insertion_detected"
#'
#' @keywords internal
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
        r_cigar  = read_stats$cigar,
        r_qual   = read_stats$qual
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
