# Project : ElsaBLab_fRagmentomics

#' Get base identity and sequencing quality from read and position of interest.
#'
#' @inheritParams get_insertion
#'
#' @return a named list with names `base` and `qual`.
#'
#' @importFrom stringr str_extract
#'
#' @noRd
get_base_qual_from_read <- function(pos, r_pos, r_cigar, r_query, r_qual) {
  # table of operations-cursor movement
  op <- c("M", "I", "D", "N", "S", "H", "P", "=", "X")
  consumes_query <- c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE)
  consumes_reference <- c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE)
  cigar_table <- data.frame(op, consumes_query, consumes_reference)

  # Cursors to track our position
  cursor_reference <- r_pos
  index_query <- 1

  while (nchar(r_cigar) > 0) {
    # Parse the CIGAR string for the next operation
    c_dg_str <- str_extract(r_cigar, "^[\\d]+")
    c_dg_num <- as.numeric(c_dg_str)
    c_op <- str_extract(r_cigar, "(?<=[\\d])[MIDNSHP=X]")
    r_cigar <- gsub(paste0("^", c_dg_str, c_op), "", r_cigar)

    # Get properties of the current operation
    op_info <- cigar_table[cigar_table$op == c_op, ]

    # Calculate how much this operation moves the cursors
    ref_move <- if (op_info$consumes_reference) c_dg_num else 0
    query_move <- if (op_info$consumes_query) c_dg_num else 0

    # Check if the target position falls within the genomic block
    # covered by THIS operation. This only applies to ops that consume the reference.
    if (op_info$consumes_reference) {
      if (pos >= cursor_reference && pos < (cursor_reference + ref_move)) {
        if (c_op %in% c("M", "=", "X")) {
          # The position is a match/mismatch. Find the exact base.
          offset <- pos - cursor_reference
          target_index <- index_query + offset
          base <- substr(r_query, target_index, target_index)
          qual <- substr(r_qual, target_index, target_index)
          return(list(base = base, qual = qual))
        } else if (c_op %in% c("D", "N")) {
          # The position is a deletion in the read relative to the reference.
          return(list(base = "-", qual = NA))
        }
      }
    }

    # If the position was not in the current block, update the cursors
    # to point to the beginning of the next block.
    cursor_reference <- cursor_reference + ref_move
    index_query <- index_query + query_move
  }

  # If the loop finishes, the read did not cover the position.
  return(list(base = NA, qual = NA))
}
