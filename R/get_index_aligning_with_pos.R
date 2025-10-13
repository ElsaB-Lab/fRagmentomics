#' Find the read sequence index for a genomic position
#'
#' @description Parses a read's CIGAR string to find the 1-based index in the
#' read's sequence that corresponds to a specific 1-based genomic coordinate.
#'
#' @inheritParams get_base_basq_mstat_from_read
#' @return An integer scalar representing the 1-based index in the read
#' sequence.
#' \itemize{
#'   \item Returns '-1' if the read does not cover the position.
#'   \item Returns '-2' if the position falls within a deletion or skipped
#'     region ('D' or 'N' CIGAR operation) in the read.
#' }
#'
#' @keywords internal
get_index_aligning_with_pos <- function(pos, read_stats) {
    # get POS, CIGAR, SEQ, and QUAL from the read
    read_pos <- read_stats$POS
    read_cigar <- read_stats$CIGAR

    # table of operations-cursor movement
    operations <- c("M", "I", "D", "N", "S", "H", "P", "=", "X")
    consumes_seq <- c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE)
    consumes_ref <- c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE)
    cigar_operations <- data.frame(operations, consumes_seq, consumes_ref)

    # Cursors to track our position along the reference and the read ref_pos
    # starts at the read's starting position (1-based) read_idx starts at the
    # first base of the query sequence (1-based)
    ref_pos <- read_pos
    read_idx <- 1

    while (nchar(read_cigar) > 0) {
        # get current cigar operation and number of associated bases
        read_cigar_nb_str <- str_extract(read_cigar, "^[\\d]+")
        read_cigar_nb <- as.numeric(read_cigar_nb_str)
        read_cigar_op <- str_extract(read_cigar, "(?<=[\\d])[MIDNSHP=X]")

        # get reduced cigar without the current operation
        read_cigar <- gsub(paste0("^", read_cigar_nb_str, read_cigar_op), "", read_cigar)

        # get properties of the current operation
        op_consumes <- cigar_operations[cigar_operations$operations == read_cigar_op,
            ]
        consumes_seq <- op_consumes$consumes_seq
        consumes_ref <- op_consumes$consumes_ref

        # Calculate how much this entire block moves the cursors
        ref_move <- if (consumes_ref)
            read_cigar_nb else 0
        seq_move <- if (consumes_seq)
            read_cigar_nb else 0

        # Check if the target position falls within the genomic block covered
        # by this operation. This only applies to ops that consume the
        # reference.
        if (consumes_ref) {
            if (pos >= ref_pos && pos < (ref_pos + ref_move)) {
                # The position of interest is within this block.
                if (read_cigar_op %in% c("M", "=", "X")) {
                    # The position is a match/mismatch. Find the exact index in
                    # the read.
                    offset <- pos - ref_pos
                    return(read_idx + offset)
                } else if (read_cigar_op %in% c("D", "N")) {
                    # The position is deleted a skipped in the read.
                    # As per original logic, return -2.
                    return(-2)
                }
            }
        }

        # The position searched is the start of a softclipped or inserted
        # region in the end of the read
        if (consumes_seq & (pos == ref_pos) & nchar(read_cigar) == 0) {
            return(-0.5)
        }

        # If the position was not in the current block, update the cursors to
        # point to the beginning of the next block.
        ref_pos <- ref_pos + ref_move
        read_idx <- read_idx + seq_move
    }

    # if the loop finishes, the read did not cover the position of interest
    return(-1)
}
