#' Extract Soft-Clipped Base Counts from CIGAR Strings
#'
#' @description This function extracts the number of soft-clipped bases at the 5' ends of 5' reads and 3' ends of 3' reads.
#' If no soft clipping is detected, the count is set to 0.
#'
#' @param cigar_5p Character string. The CIGAR string of the 5' read.
#' @param cigar_3p Character string. The CIGAR string of the 3' read.
#'
#' @return
#' A 'list' containing two integer:
#' \itemize{
#'   \item 'nb_softclip_5p': The number of soft-clipped bases at the start of
#'     the 5' read.
#'   \item 'nb_softclip_3p': The number of soft-clipped bases at the end of
#'     the 3' read.
#' }
#'
#' @keywords internal
get_fragment_bases_5p_3p_softclip <- function(cigar_5p, cigar_3p) {
    # For CIGAR in 5' : We look if the CIGAR start with a softlipping
    if (grepl("^(\\d+)S", cigar_5p)) {
        nb_softclip_5p <- as.numeric(sub("^(\\d+)S.*", "\\1", cigar_5p))
    } else {
        nb_softclip_5p <- 0
    }

    # For CIGAR in 3' : We look if the CIGAR end with a softlipping
    if (grepl("(\\d+)S$", cigar_3p)) {
        nb_softclip_3p <- as.numeric(sub(".*?(\\d+)S$", "\\1", cigar_3p))
    } else {
        nb_softclip_3p <- 0
    }

    return(list(
        nb_softclip_5p = nb_softclip_5p,
        nb_softclip_3p = nb_softclip_3p
    ))
}
