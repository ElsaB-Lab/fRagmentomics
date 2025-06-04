#' Get the informations about the presence of a insertion in the read.
#'
#' @inheritParams get_base_basq_mstat_from_read
#' @param type, choose between "INS" or "DEL"
#'
#' @return a boolean indicating TRUE if the insertion was found in the CIGAR a the position of interest, FALSE otherwise
#'
#' @keywords internal
search_for_indel_in_cigar <- function(pos, ref, alt, read_stats, type) {
  stopifnot(type %in% c("I", "D"))
  extra_info <- NULL

  # CIGAR operations
  operations <- c("M", "I", "D", "N", "S", "H", "P", "=", "X")
  consumes_seq <- c("yes", "yes", "no", "no", "yes", "no", "no", "yes", "yes")
  consumes_ref <- c("yes", "no", "yes", "yes", "no", "no", "no", "yes", "yes")
  cigar_operations <- data.frame(operations, consumes_seq, consumes_ref)

  # Parse the CIGAR string
  ops <- parse_cigar(read_stats$CIGAR)

  # Initialiaze position in reference and cursor in read sequence
  ref_pos <- read_stats$POS
  read_idx <- 0

  # Iterating over CIGAR operations
  for (i in seq_len(nrow(ops))) {
    op_len <- ops$length[i]
    op_type <- ops$type[i]

    # test if we have found the insertion we are looking for
    if (op_type==type){
      # check that we have a CIGAR operation of the right type, of the right size, at the right position
      if (op_type=="I" && ref_pos - 1 == pos){
        if (op_len==nchar(alt)-1 && substr(read_stats$SEQ, read_idx, read_idx + op_len)==alt){
          # Of note, if the CIGAR starts with I, then read_idx is zero and the substr
          # below will return a string of size op_len instead of op_len + 1
          # In that case, the sequence extracted from the read will not include the base before the insertion
          # which is included in the 'alt' sequence and the comparison below will not hold.
          return(list(TRUE, NULL))
        } else {
          return(list(FALSE, "other MUT found in CIGAR"))
        }
      } else if (op_type=="D" && ref_pos - 1 == pos){
        if (op_len==nchar(ref)-1){
          return(list(TRUE, NULL))
        } else {
          return(list(FALSE, "other MUT found in CIGAR"))
        }
      }
    }

    # execute cigar operation to move cursors
    consumes_seq <- cigar_operations[cigar_operations$op == op_type, "consumes_seq"]
    consumes_ref <- cigar_operations[cigar_operations$op == op_type, "consumes_ref"]

    # execute read_cigar operation along the reference
    if (consumes_ref == "yes") {
      ref_pos <- ref_pos + op_len
    }

    # execute read_cigar operation along the read
    if (consumes_seq == "yes") {
      read_idx <- read_idx + op_len
    }
  }

  return(list(FALSE, "MUT not in CIGAR"))
}
