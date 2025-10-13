#' Trim SoftClipped Bases from a Read
#'
#' @description This function parses a CIGAR string to identify and remove soft-clipped bases from the 5p and 3p
#' ends of a sequence and its corresponding quality string. It also updates the CIGAR string.
#'
#' @param read_stats A list containing CIGAR, SEQ, and QUAL strings.
#'
#' @return A list with the updated CIGAR, SEQ, and QUAL strings.
#'
#' @importFrom stringr str_match str_remove
#'
#' @noRd
remove_softclip <- function(read_stats) {
    cigar <- read_stats$CIGAR
    seq <- read_stats$SEQ
    qual <- read_stats$QUAL

    # Get back the number of soft clipping in 5'. Ignore a leading H.
    softclip_5p <- str_match(cigar, "^(?:\\d+H)?(\\d+)S")
    n_softclip_5p <- if (is.na(softclip_5p[1, 2]))
        0 else as.integer(softclip_5p[1, 2])

    # Get back the number of soft clipping in 3'. Ignore a leading H.
    softclip_3p <- str_match(cigar, "(\\d+)S(?:\\d+H)?$")
    n_softclip_3p <- if (is.na(softclip_3p[1, 2]))
        0 else as.integer(softclip_3p[1, 2])

    # Trim sequence and quality
    seq_len <- nchar(seq)
    new_seq <- substr(seq, start = n_softclip_5p + 1, stop = seq_len - n_softclip_3p)
    new_qual <- substr(qual, start = n_softclip_5p + 1, stop = seq_len - n_softclip_3p)

    # Update cigar
    new_cigar <- sub("^(\\d+H)?\\d+S", "\\1", cigar)
    new_cigar <- sub("\\d+S(\\d+H)?$", "\\1", new_cigar)

    # Update read length
    new_read_length <- nchar(new_seq)

    # Return the new seq, qual and cigar
    return(list(SEQ = new_seq, QUAL = new_qual, CIGAR = new_cigar, read_length = new_read_length))
}
