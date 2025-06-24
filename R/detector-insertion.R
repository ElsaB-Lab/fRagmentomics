#' Get the informations about the presence of an insertion in the read.
#'
#' @param pos numeric value representing the position of interest
#' @param alt character vector representing the insertion sequence
#' @param r_pos numeric value representing the read mapping position
#' @param r_cigar character vector representing the read CIGAR
#' @param r_query character vector representing the read base sequence
#' @param r_qual character vector representing the read sequencing qualities
#' @param pos_after_indel_repetition integer, 1-based genomic position of the
#'        first reference nucleotide after the (potentially repetitive)
#'        region involved in the indel. This position is used
#'        to determine if a read is long enough to be considered
#'        unambiguous. If the read does not cover sufficiently up to
#'        this position, it may be marked "ambiguous".
#'
#' @return a named list with names 'base','qual'.
#'         'base' can be:
#'           - 'NA': 'pos' is invalid, or read does not cover 'pos'.
#'           - '"ambiguous"': Read coverage is insufficient relative to 'pos_after_indel_repetition'
#'               to make an unambiguous call. This status is assigned if
#'               'last_ref_pos_covered_by_read < pos_after_indel_repetition'.
#'               This determination overrides CIGAR-based findings if the read
#'               is deemed too short by this rule.
#'           - '"+INSERTED_SEQ"': The specific insertion CIGAR operation was
#'             detected (e.g., '"+CAG"') and the read is long enough
#'             to be considered unambiguous.
#'           - '"no_insertion_detected"': No specific insertion CIGAR operation
#'             was found and the read is long enough to be considered
#'             unambiguous (WT).
#'         'qual' provides:
#'           - 'NA' if 'base' is 'NA', or if an insertion
#'             CIGAR operation is found at the very start of the aligned read
#'             portion (and the call is not "ambiguous").
#'           - "ambiguous" if base is "ambiguous".
#'           - The quality score of the read base preceding the insertion.
#'           - "no_insertion_detected" if base is "no_insertion_detected".
#'
#' @keywords internal
get_insertion <- function(pos, alt, r_pos, r_cigar, r_query, r_qual, pos_after_indel_repetition) {
  # Initialize variables
  c_base_found <- NULL
  c_qual_found <- NULL
  insertion_op_found_in_cigar <- FALSE

  # Length of the purely inserted sequence
  insertion_actual_length <- nchar(alt) - 1

  if (insertion_actual_length <= 0) {
    warning(paste0(
      "'get_insertion': 'alt' sequence (insertion) has length ", insertion_actual_length,
      ". Expected a positive length for the inserted sequence. No insertion will be detected."
    ))
  }

  ops <- parse_cigar(r_cigar)
  read_cursor <- 0 # 0-based count of consumed read bases
  current_ref_genome_pos <- r_pos # 1-based current position in the reference

  for (i in seq_len(nrow(ops))) {
    op_len <- ops$length[i]
    op_type <- ops$type[i]

    ref_pos_at_start_of_op <- current_ref_genome_pos

    if (op_type %in% c("M", "=", "X")) {
      current_ref_genome_pos <- current_ref_genome_pos + op_len
      read_cursor <- read_cursor + op_len
    } else if (op_type %in% c("N", "D")) {
      current_ref_genome_pos <- current_ref_genome_pos + op_len
    } else if (op_type == "I") {
      # Check if this is the target insertion operation
      # 'pos' is the 1-based ref coord of base before insertion.
      # 'current_ref_genome_pos' is the ref coord at which the insertion occurs (i.e., after 'pos').
      # So, current_ref_genome_pos - 1 should equal 'pos'.
      if ((ref_pos_at_start_of_op - 1 == pos) &&
        op_len == insertion_actual_length) {
        # Extract inserted sequence from r_query.
        if (read_cursor > 0 &&
          (substr(r_query, read_cursor, read_cursor) == substr(alt, 1, 1))) {
          seq_from_read <- substr(r_query, read_cursor, read_cursor + op_len)

          insertion_op_found_in_cigar <- TRUE
          c_base_found <- paste0("+", substring(seq_from_read, 2))

          # Extract the quality of the inserted sequence
          c_qual_found <- substr(r_qual, read_cursor, read_cursor)
        } else if (read_cursor > 0 &&
          (substr(r_query, read_cursor, read_cursor) != substr(alt, 1, 1))) {
          insertion_op_found_in_cigar <- TRUE
          c_base_found <- "no_insertion_detected"
          c_qual_found <- "no_insertion_detected"
        } else {
          insertion_op_found_in_cigar <- TRUE
          c_base_found <- NA_character_ # Insertion at the  start of the read
          c_qual_found <- NA_character_
        }
      }
      read_cursor <- read_cursor + op_len # 'I' op consumes read bases
    } else if (op_type == "S") {
      read_cursor <- read_cursor + op_len # 'S' op consumes read bases (soft clip)
    } else if (op_type %in% c("H", "P")) {
      # H (Hard clip), P (Padding) - No change to read_cursor for aligned part, no change to ref pos.
    }
  } # End of CIGAR processing loop

  last_ref_pos_covered_by_read <- current_ref_genome_pos - 1

  #-------------------------------------------------
  # Final decision on insertion status
  #-------------------------------------------------
  # Handle invalid input.
  if (pos < 1) {
    warning(paste0(
      "'get_insertion': 'pos' (", pos,
      ") is not a valid 1-based coordinate for the base preceding an insertion. Returning NA."
    ))
    return(list(base = NA_character_, qual = NA_character_))
  }

  # Check for ambiguity, as this overrides all other findings.
  can_evaluate_ambiguity_rule <- TRUE
  if (is.null(pos_after_indel_repetition) || is.na(pos_after_indel_repetition) || !is.numeric(pos_after_indel_repetition)) {
    warning("'get_insertion': 'pos_after_indel_repetition' is missing, NA, or not numeric. Ambiguity rule based on read length cannot be applied.")
    can_evaluate_ambiguity_rule <- FALSE
  } else if (pos_after_indel_repetition <= 0) {
    warning(paste0(
      "'get_insertion': 'pos_after_indel_repetition' (", pos_after_indel_repetition,
      ") must be positive for the ambiguity rule. Rule cannot be applied meaningfully."
    ))
    can_evaluate_ambiguity_rule <- FALSE
  }

  # Ambiguity rule is also not meaningful if there's no valid insertion length.
  if (can_evaluate_ambiguity_rule && insertion_actual_length <= 0) {
    can_evaluate_ambiguity_rule <- FALSE
  }

  # If the read is NOT ambiguous, then check if it covers the insertion site.
  if (last_ref_pos_covered_by_read < pos || r_pos > pos) {
    # The read does not cover the position of interest.
    return(list(base = NA_character_, qual = NA_character_))
  }

  if (can_evaluate_ambiguity_rule) {
    if (last_ref_pos_covered_by_read < pos_after_indel_repetition) {
      # The read is too short to make a definitive call, so it is ambiguous.
      return(list(base = "ambiguous", qual = "ambiguous"))
    }
  }

  # If the read is long enough and covers the site, return the call
  # based on the CIGAR string analysis performed earlier.
  if (insertion_op_found_in_cigar) {
    return(list(base = c_base_found, qual = c_qual_found))
  } else {
    return(list(base = "no_insertion_detected", qual = "no_insertion_detected"))
  }
}
