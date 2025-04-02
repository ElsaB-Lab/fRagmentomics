# Project : ElsaBLab_fRagmentomics

#' Retrieve mutation-specific information:
#' Extracts base, quality, and indel informations based on the mutation type.
#'
#' @inheritParams process_fragment
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
        r_pos    = read_stats$POS,
        r_cigar  = read_stats$CIGAR,
        r_query  = read_stats$SEQ,
        r_qual   = read_stats$QUAL
      )
    )
  } else if (mutation_type == "deletion") {
    return(
      get_deletion(
        pos      = pos,
        ref      = ref,
        r_pos    = read_stats$POS,
        r_cigar  = read_stats$CIGAR,
        r_qual   = read_stats$QUAL
      )
    )
  } else if (mutation_type == "insertion") {
    return(
      get_insertion(
        pos      = pos,
        alt      = alt,
        r_pos    = read_stats$POS,
        r_cigar  = read_stats$CIGAR,
        r_query  = read_stats$SEQ,
        r_qual   = read_stats$QUAL
      )
    )
  } else {
    return(list(
      base = "Error: mutation_type not defined properly",
      qual = "Error: mutation_type not defined properly"
    ))
  }
}
