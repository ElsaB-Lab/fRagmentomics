# Project : ElsaBLab_fRagmentomics

#' Extract Fragment Bases and Quality Scores from Read Sequences
#' This function extracts a fixed number of bases (`nbr_bases`)
#' and their corresponding quality scores from the 5' and 3' ends of reads.
#'
#' @param nbr_bases The number of bases to extract from each end (5' and 3').
#' @param query_5p The nucleotide sequence of the 5' read.
#' @param query_3p The nucleotide sequence of the 3' read.
#' @param qual_5p The quality string corresponding to the 5' read.
#' @param qual_3p The quality string corresponding to the 3' read.
#'
#' @return A list containing four character elements:
#'   -fragment_bases_5p: The first `nbr_bases` nucleotides from the 5' read.
#'   -fragment_bases_3p: The last `nbr_bases` nucleotides from the 3' read.
#'   -fragment_Qbases_5p: The first `nbr_bases` qual chr from the 5' read.
#'   -fragment_Qbases_3p: The last `nbr_bases` qual chr from the 3' read.
#'
#' @keywords internal
process_get_fragment_bases <- function(
    nbr_bases,
    query_5p,
    query_3p,
    qual_5p,
    qual_3p) {
  if (nbr_bases > nchar(query_5p)) {
    warning("Number of bases in 5p and 3p bigger than the fragment length")
    fragment_bases_5p <- query_5p
    fragment_qbases_5p <- qual_5p
    fragment_bases_3p <- query_3p
    fragment_qbases_3p <- qual_3p
  } else {
    fragment_bases_5p <- substr(query_5p, 1, nbr_bases)
    fragment_qbases_5p <- substr(qual_5p, 1, nbr_bases)
    fragment_bases_3p <- substr(
      query_3p, nchar(query_3p) - nbr_bases + 1, nchar(query_3p)
    )
    fragment_qbases_3p <- substr(
      qual_3p, nchar(qual_3p) - nbr_bases + 1, nchar(qual_3p)
    )
  }

  return(list(
    fragment_bases_5p = fragment_bases_5p,
    fragment_bases_3p = fragment_bases_3p,
    fragment_Qbases_5p = fragment_qbases_5p,
    fragment_Qbases_3p = fragment_qbases_3p
  ))
}
