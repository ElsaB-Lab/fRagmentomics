# Project : ElsaBLab_fRagmentomics

#' Get the informations about the presence of a deletion in the read.
#'
#' @inheritParams get_insertion
#' @param ref character vector representing the VCF REF allele for the deletion.
#' @param pos_after_indel_repetition integer, 1-based genomic position of the
#' first reference nucleotide after the (potentially repetitive) region involved
#' in the indel. This position is used to determine if a read not showing the
#' deletion is long enough to be called Wild Type or if it's ambiguous.
#'
#' @return a named list with names 'base','qual'.
#' 'base' can be:
#'    - "no_deletion_detected": Read covers the site, shows no deletion, and is
#'      long enough to be unambiguous = WT.
#'    - "-DELETED_SEQ": The specific deletion was detected (e.g., "-ATC").
#'    - NA_character : Read does not cover the deletion site.
#'    - "ambiguous": Read coverage is insufficient relative to 'pos_after_indel_repetition'
#'      to make an unambiguous call. This status is assigned if
#'      'last_ref_pos_covered_by_read < (pos_after_indel_repetition - deletion_length)'.
#' 'qual' provides the quality score of the base preceding the deletion if
#'    detected, "ambiguous" if base is "ambiguous", NA if base is
#'    NA, or "no_deletion_detected".
#'
#' @keywords internal
get_info_deletion <- function(pos, ref, r_pos, r_cigar, r_qual, pos_after_indel_repetition) {
  # Initialize variables for storing if the specific deletion CIGAR op is found
  c_base_found <- NULL # Sequence like "-ATC" if D op found
  c_qual_found <- NULL # Associated quality if D op found
  deletion_op_found_in_cigar <- FALSE # Flag: TRUE if the specific D op is found

  # Calculate the length of the deleted sequence based on VCF REF convention
  deletion_length <- nchar(ref) - 1

  # Parse the CIGAR string
  ops <- parse_cigar(r_cigar)
  read_cursor <- 0 # 0-based cursor for current position in the read
  current_ref_genome_pos <- r_pos # 1-based current position in the reference

  # Iterate over CIGAR operations to:
  # 1. Check for the specific deletion CIGAR operation.
  # 2. Determine the last position covered by the read ('last_ref_pos_covered_by_read').
  for (i in seq_len(nrow(ops))) {
    op_len <- ops$length[i]
    op_type <- ops$type[i]

    # ref_pos_at_start_of_op is the 1-based reference genome position
    # where the current CIGAR operation begins.
    ref_pos_at_start_of_op <- current_ref_genome_pos

    if (op_type %in% c("M", "=", "X")) {
      current_ref_genome_pos <- current_ref_genome_pos + op_len
      read_cursor <- read_cursor + op_len
    } else if (op_type %in% c("N", "D")) {
      # Check if this is the target deletion operation
      if (op_type == "D" && (ref_pos_at_start_of_op - 1 == pos) && op_len == deletion_length) {
        deletion_op_found_in_cigar <- TRUE
        c_base_found <- paste0("-", substring(ref, 2))

        if (read_cursor > 0) { # read_cursor is count of preceding read bases
          # substr is 1-based in R. If read_cursor = 3, it means 3 bases came before.
          # So, quality of 3rd base (last matched base) is substr(r_qual, read_cursor, read_cursor).
          c_qual_found <- substr(r_qual, read_cursor, read_cursor)
        } else {
          # No preceding base in the read to take quality from.
          c_qual_found <- NA_character_
        }
      }
      current_ref_genome_pos <- current_ref_genome_pos + op_len # Advance ref pos for D or N
    } else if (op_type %in% c("I", "S")) {
      read_cursor <- read_cursor + op_len # Consumes read bases
    } else if (op_type %in% c("H", "P")) {
      # H (Hard clip), P (Padding) - No change to read_cursor for aligned part, no change to ref pos.
    }
  } # End of CIGAR processing loop

  # 'current_ref_genome_pos' is now the 1-based ref coord after the last ref-consuming op.
  # 'last_ref_pos_covered_by_read' is the last ref base (1-based) covered by alignment.
  last_ref_pos_covered_by_read <- current_ref_genome_pos - 1

  #-------------------------------------------------
  # Final decision on deletion status
  #-------------------------------------------------

  # Check if the position is valid.
  if (pos < 1) {
    warning(paste0(
      "'get_deletion': 'pos' (", pos,
      ") is not a valid 1-based coordinate for the base preceding a deletion. Returning NA."
    ))
    return(list(base = NA_character_, qual = NA_character_))
  }

  # Check if the read covers the indel position.
  if (last_ref_pos_covered_by_read < pos) {
    # The read does not cover the position of interest.
    return(list(base = NA_character_, qual = NA_character_))
  }

  # Check for ambiguity.
  can_evaluate_ambiguity_rule <- TRUE
  if (is.null(pos_after_indel_repetition) || is.na(pos_after_indel_repetition) || !is.numeric(pos_after_indel_repetition)) {
    warning("'get_deletion': 'pos_after_indel_repetition' is missing, NA, or not numeric. Ambiguity rule based on read length cannot be applied.")
    can_evaluate_ambiguity_rule <- FALSE
  } else if (pos_after_indel_repetition <= 0) {
    warning(paste0(
      "'get_deletion': 'pos_after_indel_repetition' (", pos_after_indel_repetition,
      ") must be positive for the ambiguity rule. Rule cannot be applied meaningfully."
    ))
    can_evaluate_ambiguity_rule <- FALSE
  }

  # The ambiguity rule cannot be applied if it's not a positive-length deletion.
  if (can_evaluate_ambiguity_rule && deletion_length <= 0) {
    warning(paste0(
      "'get_deletion': 'deletion_length' (", deletion_length,
      ") must be positive for the ambiguity rule related to deletions. Rule cannot be applied in this context."
    ))
    can_evaluate_ambiguity_rule <- FALSE
  }

  if (can_evaluate_ambiguity_rule) {
    threshold_ref_pos <- pos_after_indel_repetition - deletion_length

    if (nrow(ops) > 0) { # Ensure the CIGAR string was not empty
      if (op_type == "D") {
        last_ref_pos_covered_by_read <- last_ref_pos_covered_by_read - op_len
      }
    }

    if (last_ref_pos_covered_by_read < threshold_ref_pos) {
      # The read is too short to be certain, so it's ambiguous. Return immediately.
      return(list(base = "ambiguous", qual = "ambiguous"))
    }
  }

  # If the read is long enough (not ambiguous) and covers the site, return the result
  # based on the CIGAR string analysis performed earlier.
  if (deletion_op_found_in_cigar) {
    return(list(base = c_base_found, qual = c_qual_found))
  } else {
    # No D op found, not NA, not ambiguous by length -> "no_deletion_detected".
    return(list(base = "no_deletion_detected", qual = "no_deletion_detected"))
  }
}
