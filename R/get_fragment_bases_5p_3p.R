#' Extract Fragment Bases and Quality Scores from Read Sequences
#'
#' @description This function extracts a fixed number of bases ('n_bases')
#' and their corresponding quality scores from the 5' and 3' ends of reads.
#'
#' @param n_bases The number of bases to extract from each end (5' and 3').
#' @param seq_5p The nucleotide sequence of the 5' read.
#' @param seq_3p The nucleotide sequence of the 3' read.
#' @param qual_5p The quality string corresponding to the 5' read.
#' @param qual_3p The quality string corresponding to the 3' read.
#'
#' @return A list containing four character elements:
#'   -fragment_bases_5p: The first 'n_bases' nucleotides from the 5' read.
#'   -fragment_bases_3p: The last 'n_bases' nucleotides from the 3' read.
#'   -fragment_basqs_5p: The first 'n_bases' qual chr from the 5' read.
#'   -fragment_basqs_3p: The last 'n_bases' qual chr from the 3' read.
#'
#' @keywords internal
get_fragment_bases_5p_3p <- function(n_bases, seq_5p, seq_3p, qual_5p, qual_3p) {
  if (n_bases > nchar(seq_5p)) {
    warning("Number of bases in 5p and 3p bigger than the fragment length")
    fragment_bases_5p <- seq_5p
    fragment_basqs_5p <- qual_5p
    fragment_bases_3p <- seq_3p
    fragment_basqs_3p <- qual_3p
  } else {
    fragment_bases_5p <- substr(seq_5p, 1, n_bases)
    fragment_basqs_5p <- substr(qual_5p, 1, n_bases)
    fragment_bases_3p <- substr(
      seq_3p, nchar(seq_3p) - n_bases + 1,
      nchar(seq_3p)
    )
    fragment_basqs_3p <- substr(
      qual_3p, nchar(qual_3p) - n_bases + 1,
      nchar(qual_3p)
    )
  }

  return(list(
    fragment_bases_5p = fragment_bases_5p,
    fragment_bases_3p = fragment_bases_3p,
    fragment_basqs_5p = fragment_basqs_5p,
    fragment_basqs_3p = fragment_basqs_3p
  ))
}
