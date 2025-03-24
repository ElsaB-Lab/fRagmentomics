# Project : ElsaBLab_fRagmentomics

#' Get the informations about the presence of a deletion in the read.
#'
#' @inheritParams fRagmentomics
#' @param pos numeric value representing the position of interest
#' @param ref character vector representing the deletion sequence
#' @param r_pos numeric value representing the read mapping position
#' @param r_cigar character vector representing the read CIGAR
#' @return a named list with names `base`,`qual`
#'
#' @noRd
get_deletion <- function(pos, ref, r_pos, r_cigar) {
  # Initialize variables
  c_base <- "no_deletion_detected"
  c_qual <- "no_deletion_detected"                 
  ref_len <- nchar(ref) -1      # because we added one nucleotide for VCF convention
  pos_is_readed <- FALSE

  # Parse the CIGAR string
  # ops is a df with a column length and a column type; Each row for each operation
  ops <- parse_cigar(r_cigar)

  # Iterating over CIGAR operations
  current_pos <- r_pos  # Current reference position, aligned with the read
  for (i in seq_len(nrow(ops))) {
    op_len <- ops$length[i]
    op_type <- ops$type[i]

    if (op_type %in% c("M", "N", "=", "X")) {
      # M: match or mismatch, N: skip in the reference, =: perfect match, X: mismatch
      current_pos <- current_pos + op_len

    } else if (op_type == "D") {
      # D: deletion in the read (missing in the read compared to the reference)
      # Compare this deletion to the one in the reference
      diff_mutation_bam <- pos - (current_pos - 1)
      
      # Advance the reference position after the deletion
      current_pos <- current_pos + op_len

      # Check if position and length of deletion is correct
      # Deletion sequence has already been check in the preprocessing
      if ((ref_len == op_len) && (diff_mutation_bam == 0)) {
        c_base <- "deletion_detected"
        c_qual <- "-"
        pos_is_readed <- TRUE
        break
      }

    } else {
      # I, S, H, P: Insertion, Soft clip, Hard clip, Pad - do not advance the reference position
      # No action is taken
    }
  }

  # Check if read covers the position of interest
  if ((current_pos - 1) > pos) {
    pos_is_readed <- TRUE
  }

  # if the read does not cover the position of interest
  if (pos_is_readed == FALSE){
    c_base <- NA
    c_qual <- NA
  }

  return(list(base = c_base, qual = c_qual))
}