# Project : ElsaBLab_fRagmentomics

#' Get the informations about the presence of a insertion in the read.
#'
#' @inheritParams build_fragments_info_table*
#' @param pos numeric value representing the position of interest
#' @param alt character vector representing the insertion sequence
#' @param r_pos numeric value representing the read mapping position
#' @param r_cigar character vector representing the read CIGAR
#' @param r_query character vector representing the read base sequence
#' @param r_qual character vector representing the read sequencing qualities
#' @return a named list with names `base`,`qual` and `indel`.
#'
#' @importFrom stringr str_extract utils
#'
#' @noRd
get_insertion <- function(pos, alt, r_pos, r_cigar, r_query, r_qual) {
  c_base <- NA
  c_qual <- NA
  c_indel <- "No_support_ins"  # By default, no insertion
  alt_len <- nchar(alt)

  # Parse the CIGAR string
  ops <- parse_cigar(r_cigar)

  # Iterating over CIGAR operations
  current_pos <- r_pos  # Current reference position, aligned with the read
  for (i in seq_len(nrow(ops))) {
    op_len <- ops$length[i]
    op_type <- ops$type[i]

    if (op_type == "M") {
      # M: match or mismatch (but advances the reference position)
      current_pos <- current_pos + op_len
    } else if (op_type == "D") {
      # D: deletion in the read (advances the reference position but not the read)
      # Therefore, we modify current_pos
      current_pos <- current_pos + op_len
    } else if (op_type == "I") {
      # I: insertion in the read (addition in the read compared to the reference)
      # Compare this insertion with the reference
      diff_mutation_bam <- pos - (current_pos - 1)
      is_length_ok <- (alt_len == op_len)

      ins_rep <- define_ins_rep(alt_len, current_pos, r_pos, r_cigar, r_query)

      c_indel <- check_seq_rep(ins_rep, is_length_ok, diff_mutation_bam, op_len)

      if (c_indel == 1) {
        c_qual <- find_c_qual(alt_len, current_pos, r_pos, r_cigar, r_qual)
        # A corresponding insertion was found
        break
      }
    } else if (op_type %in% c("N", "=", "X")) {
      # N: skip in the reference, =: perfect match, X: mismatch
      # All advance the reference position
      current_pos <- current_pos + op_len
    } else {
      # S, H, P: Soft clip, Hard clip, Pad - do not advance the reference position
      # No action is taken
    }
  }

  return(list(base = c_base, qual = c_qual, indel = c_indel))
}