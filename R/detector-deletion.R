# Project : ElsaBLab_fRagmentomics

#' Get the informations about the presence of a deletion in the read.
#'
#' @inheritParams build_fragments_info_table*
#' @param chr numeric value representing the chromosome
#' @param pos numeric value representing the position of interest
#' @param ref character vector representing the deletion sequence
#' @param del_info character vector sous la forme "$pos_final,$rep_del"
#' @param r_pos numeric value representing the read mapping position
#' @param r_cigar character vector representing the read CIGAR
#' @return a named list with names `base`,`qual` and `indel`.
#'
#' @noRd
get_deletion <- function(chr, pos, ref, del_info, r_pos, r_cigar) {
  # Split the del_info parameter into final_pos and del_rep
  vars <- strsplit(del_info, ",")[[1]]
  final_pos <- suppressWarnings(as.numeric(vars[1]))
  del_rep <- as.numeric(vars[2])

  # Initialize variables
  c_base <- NA
  c_qual <- NA
  c_indel <- "No_support_del"   # by default, no deletion
  ref_len <- nchar(ref)

  # Check if final_pos is not numeric
  if (is.na(final_pos)) {
    c_indel <- "no_del_found_in_ref_genome"
    return(list(base = c_base, qual = c_qual, indel = c_indel))

  } else {
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
      } else if (op_type == "I") {
        # I: insertion in the read (does not advance the reference position, only advances in the read)
        # Therefore, we do not modify current_pos
        # (current_pos remains the reference position)
      } else if (op_type == "D") {
        # D: deletion in the read (missing in the read compared to the reference)
        # Compare this deletion to the one in the reference
        diff_mutation_bam <- final_pos - (current_pos - 1)
        is_length_ok <- (ref_len == op_len)

        c_indel <- check_seq_rep(del_rep, is_length_ok, diff_mutation_bam, op_len)

        # Advance the reference position after the deletion
        current_pos <- current_pos + op_len

        if (c_indel == 1) {
          # A corresponding deletion was found
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
}