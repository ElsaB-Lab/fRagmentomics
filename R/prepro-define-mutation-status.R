# Project: ElsaBLab_fRagmentomics

#' Define Mutation Status
#'
#' Determines the type of mutation (SNV, MNP, insertion, or deletion) based on
#' the length of the REF and ALT strings.
#'
#' @param ref Reference allele (character string).
#' @param alt Alternate allele (character string).
#'
#' @return A character string indicating the type of mutation: 
#'   one of "SNV", "MNP", "insertion", or "deletion".
#'
#' @noRd
define_mutation_status <- function(ref, alt) {
    # Determine the type of mutation
    mutation_status <- ""

    if (nchar(ref) == nchar(alt)) {
        if (nchar(ref) == 1) {
            mutation_status <- "SNV"     
        } else {
            mutation_status <- "MNP"       
        }
    } else if (nchar(ref) > nchar(alt)) {
        mutation_status <- "deletion"  
    } else {
        mutation_status <- "insertion"  
    }

    return(mutation_status)
}
