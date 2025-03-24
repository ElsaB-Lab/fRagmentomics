# Project : ElsaBLab_fRagmentomics

#' Extract Soft-Clipped Base Counts from CIGAR Strings
#'
#' This function extracts the number of soft-clipped bases at the 5' and 3' ends of sequencing reads
#' based on their CIGAR strings. If soft clipping is reported, it checks whether the CIGAR string starts
#' (5' end) or ends (3' end) with a soft-clipped segment (e.g., "45S"), and extracts the numeric value.
#' If no soft clipping is detected, the count is set to 0.
#'
#' @param cigar_5p Character string. The CIGAR string of the 5' read.
#' @param cigar_3p Character string. The CIGAR string of the 3' read.
#'
#' @return A list with two numeric elements:
#'  -Nb_SC_5p: Number of soft-clipped bases at the 5' end.
#'  -Nb_SC_3p: Number of soft-clipped bases at the 3' end.
#'
#' @noRd
process_get_fragment_bases_softclip <- function(cigar_5p, cigar_3p){
    # For CIGAR in 5' : We look if the CIGAR start with a softlipping
    if (grepl("^(\\d+)S", cigar_5p)) {
        Nb_SC_5p <- as.numeric(sub("^(\\d+)S.*", "\\1", cigar_5p))
    } else {
        Nb_SC_5p <- 0
    }
      
    # For CIGAR in 3' : We look if the CIGAR end with a softlipping
    if (grepl("(\\d+)S$", cigar_3p)) {
        Nb_SC_3p <- as.numeric(sub(".*?(\\d+)S$", "\\1", cigar_3p))
    } else {
        Nb_SC_3p <- 0
    }
  
    return(list(Nb_SC_5p = Nb_SC_5p, 
                Nb_SC_3p = Nb_SC_3p))
}