# Project : ElsaBLab_fRagmentomics

#' Extract Fragment Bases and Quality Scores from Read Sequences
#'
#' This function extracts a fixed number of bases (`nbr_bases`) and their corresponding quality scores 
#' from the 5' and 3' ends of sequencing reads. It takes the base sequences and quality strings from 
#' both ends of a paired-end fragment and returns the specified number of nucleotides and quality scores.
#'
#' @param nbr_bases Integer. The number of bases to extract from each end (5' and 3').
#' @param query_5p Character string. The nucleotide sequence of the 5' read.
#' @param query_3p Character string. The nucleotide sequence of the 3' read.
#' @param qual_5p Character string. The quality string corresponding to the 5' read.
#' @param qual_3p Character string. The quality string corresponding to the 3' read.
#'
#' @return A list containing four character elements:
#'   -Fragment_bases_5p: The first `nbr_bases` nucleotides from the 5' read.
#'   -Fragment_bases_3p: The last `nbr_bases` nucleotides from the 3' read.
#'   -Fragment_Qbases_5p: The first `nbr_bases` quality characters from the 5' read.
#'   -Fragment_Qbases_3p: The last `nbr_bases` quality characters from the 3' read.
#'
#' @noRd
process_get_fragment_bases <- function(nbr_bases, query_5p, query_3p, qual_5p, qual_3p){
    if (nbr_bases > nchar(query_5p)) {
        Fragment_bases_5p  <- "Number of bases in 5p and 3p bigger than the fragment length"
        Fragment_Qbases_5p <- "Number of bases in 5p and 3p bigger than the fragment length"
        Fragment_bases_3p  <- "Number of bases in 5p and 3p bigger than the fragment length"
        Fragment_Qbases_3p <- "Number of bases in 5p and 3p bigger than the fragment length"
    } else {
        Fragment_bases_5p  <- substr(query_5p, 1, nbr_bases)
        Fragment_Qbases_5p <- substr(qual_5p, 1, nbr_bases)
        Fragment_bases_3p  <- substr(query_3p, nchar(query_3p) - nbr_bases + 1, nchar(query_3p))
        Fragment_Qbases_3p <- substr(qual_3p,  nchar(qual_3p)  - nbr_bases + 1, nchar(qual_3p))
    }

    return(list(Fragment_bases_5p = Fragment_bases_5p, 
                Fragment_bases_3p = Fragment_bases_3p, 
                Fragment_Qbases_5p = Fragment_Qbases_5p, 
                Fragment_Qbases_3p = Fragment_Qbases_3p))
}