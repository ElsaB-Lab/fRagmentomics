#' Parse a CIGAR string into its component operations
#'
#' @description
#' This internal helper function deconstructs a CIGAR string into a structured
#' format, separating the operation type (e.g., M, I, D) from its length.
#'
#' @param cigar a character vector
#'
#' @return a df with length and type of the CIGAR string. One row per operation.
#'
#' @noRd
parse_cigar <- function(cigar) {
    # Extract the operations and lengths from the CIGAR string
    matches <- gregexpr("([0-9]+)([MIDNSHP=X])", cigar, perl = TRUE)
    ops_str <- regmatches(cigar, matches)[[1]]

    lengths <- as.numeric(gsub("([0-9]+)([MIDNSHP=X])", "\\1", ops_str))
    types <- gsub("([0-9]+)([MIDNSHP=X])", "\\2", ops_str)

    data.frame(length = lengths, type = types, stringsAsFactors = FALSE)
}

#' Calculate read length excluding trailing soft clips
#'
#' @description
#' This function computes the effective length of a read by subtracting the
#' number of bases that are soft-clipped at its 3' end, which is a common
#' step in variant analysis to avoid alignment artifacts.
#'
#' @param cigar A CIGAR string
#' @param seq The corresponding DNA sequence
#'
#' @return The length of the sequence after subtracting the number of bases
#'         that are soft-clipped at the end of the read.
#'
#' @noRd
calculate_len_without_end_softclip <- function(cigar, seq) {
    read_len <- nchar(seq)
    softclip_len <- 0

    # Check if the CIGAR string ends with a soft clip and capture the number of
    # bases
    if (grepl("\\d+S$", cigar)) {
        # Extract the numeric value of the trailing soft clip
        softclip_len <- as.numeric(sub(".*?(\\d+)S$", "\\1", cigar))
    }

    # Return the total length minus the trailing soft clip length
    return(read_len - softclip_len)
}

#' Extract insertion and deletion positions from a read
#'
#' @description
#' Parses a read's CIGAR string and starting position to determine the precise
#' genomic coordinates of any insertions or deletions.
#'
#' @param read_stats A list containing read infos.
#' @return A list containing two vectors: deletions and insertions.
#'
#' @noRd
get_pos_indels_from_read <- function(read_stats) {
    list_pos_del <- c()
    list_pos_ins <- c()

    # Parse the CIGAR string
    cigar_ops <- parse_cigar(read_stats$CIGAR)

    # Initialize the current genomic position
    current_ref_pos <- read_stats$POS

    # Iterate over each CIGAR operation
    for (i in seq_len(nrow(cigar_ops))) {
        op_len <- cigar_ops$length[i]
        op_type <- cigar_ops$type[i]

        # Handle operations that consume the reference
        if (op_type %in% c("M", "=", "X")) {
            # Match or Mismatch: advances the position on the reference
            current_ref_pos <- current_ref_pos + op_len
        } else if (op_type == "D") {
            # Deletion: does not advance on the read, but does on the reference
            # The positions of the deleted bases
            deleted_positions <- seq(from = current_ref_pos, length.out = op_len)
            list_pos_del <- c(list_pos_del, deleted_positions)

            # Update the position on the reference
            current_ref_pos <- current_ref_pos + op_len
        } else if (op_type == "I") {
            # Insertion: advances on the read, but not on the reference The
            # genomic position is the one just before the insertion
            pos_before_insertion <- current_ref_pos - 1

            # Create the custom positions for each inserted base
            scale <- 10^(ceiling(log10(op_len + 1)))
            inserted_positions <- pos_before_insertion + (seq_len(op_len) / scale)
            list_pos_ins <- c(list_pos_ins, inserted_positions)
            # current_ref_pos does not change
        } else if (op_type == "N") {
            # Skipped region: advances the position on the reference
            current_ref_pos <- current_ref_pos + op_len
        }
        # 'S' (soft clip), 'H' (hard clip), and 'P' (padding) operations do not
        # consume the reference, so we do nothing.
    }

    return(list(deletions = list_pos_del, insertions = list_pos_ins))
}
