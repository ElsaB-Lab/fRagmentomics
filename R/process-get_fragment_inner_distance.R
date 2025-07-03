#' Compute fragment inner distance
#' Calculates the inner distance between two paired-end sequencing reads.
#'
#' @param r_pos_5p numeric value representing the read 5p mapping position
#' @param r_pos_3p numeric value representing the read 3p mapping position
#' @param r_cigar_5p character vector representing the read 5p CIGAR
#' @param r_cigar_3p character vector representing the read 3p CIGAR
#'
#' @return a named list with names `a` and `b`.
#'
#' @importFrom stringr str_extract str_extract_all
#'
#' @keywords internal
get_fragment_inner_distance <- function(
    r_pos_5p,
    r_cigar_5p,
    r_pos_3p,
    r_cigar_3p) {
  # bases matched read 1/2
  bases_match_5p <- sum(as.numeric(
    stringr::str_extract_all(r_cigar_5p, "[[:digit:]]+(?=M)", simplify = TRUE)
  ))

  # bases deleted read 1/2
  bases_del_5p <- sum(as.numeric(
    stringr::str_extract_all(r_cigar_5p, "[[:digit:]]+(?=D)", simplify = TRUE)
  ))

  # bases inserted read 1/2
  bases_ins_5p <- sum(as.numeric(
    stringr::str_extract_all(r_cigar_5p, "[[:digit:]]+(?=I)", simplify = TRUE)
  ))

  # bases soft-clipped left read 1/2
  bases_softcl_left_3p <- as.numeric(str_extract(r_cigar_3p, "^[[:digit:]]+(?=S)"))
  bases_softcl_left_3p <- ifelse(is.na(bases_softcl_left_3p), 0, bases_softcl_left_3p)

  # bases soft-clipped right read 1/2
  bases_softcl_right_5p <- as.numeric(str_extract(r_cigar_5p, "[[:digit:]]+(?=S$)"))
  bases_softcl_right_5p <- ifelse(
    is.na(bases_softcl_right_5p), 0, bases_softcl_right_5p
  )

  # inner distance
  inner_distance <- (r_pos_3p - bases_softcl_left_3p) - (r_pos_5p + bases_match_5p + bases_del_5p + bases_softcl_right_5p - 1) - 1

  return(inner_distance)
}
