# Project : ElsaBLab_fRagmentomics

#' Get the informations about the presence of a deletion in the read.
#'
#' @inheritParams get_insertion
#' @param ref character vector representing the deletion sequence
#'
#' @return a named list with names `base`,`qual`.
#'
#' @keywords internal
get_deletion <- function(pos, ref, r_pos, r_cigar, r_qual) {
  # Initialize variables
  c_base <- "no_deletion_detected"
  c_qual <- "no_deletion_detected"
  ref_len <- nchar(ref) - 1 # because we added one nucleotide for VCF convention
  pos_is_readed <- FALSE

  # Parse the CIGAR string
  # ops is a df with columns length and type; Each row for one operation
  ops <- parse_cigar(r_cigar)
  # Current position in the lecture
  read_cursor <- 0

  # Iterating over CIGAR operations
  current_pos <- r_pos # Current reference position, aligned with the read
  for (i in seq_len(nrow(ops))) {
    op_len <- ops$length[i]
    op_type <- ops$type[i]

    if (op_type %in% c("M", "=", "X")) {
      # M: match or mismatch, N: skip in the ref, =: perfect match, X: mismatch
      current_pos <- current_pos + op_len
      read_cursor <- read_cursor + op_len
    } else if (op_type %in% c("N", "D")) {
      # D: deletion in the read (missing in the read compared to the reference)
      # Compare this deletion to the one in the reference
      if (op_type == "D" && current_pos - 1 == pos && op_len == ref_len) {
        c_base <- paste0("-", substring(ref, 2))
        c_qual <- substr(r_qual, read_cursor, read_cursor)
        pos_is_readed <- TRUE
        break
      }
      # Advance the reference position after the deletion
      current_pos <- current_pos + op_len
    } else if (op_type %in% c("I", "S")) {
      # I : insertion ; S : soft-clip  -> Consume the lecture
      read_cursor <- read_cursor + op_len
    } else {
      # H, P: Insertion, Soft clip, Hard clip, Pad
      # No advance the reference position
    }
  }

  # if the read does not cover the position of interest
  if (!pos_is_readed && (current_pos - 1) < pos) {
    c_base <- NA
    c_qual <- NA
  }

  list(base = c_base, qual = c_qual)
}
