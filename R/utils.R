#' Function to parse the CIGAR string
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

#' Calculates the length of a read without the trailing soft clipping.
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

  # Check if the CIGAR string ends with a soft clip and capture the number of bases
  if (grepl("\\d+S$", cigar)) {
    # Extract the numeric value of the trailing soft clip
    softclip_len <- as.numeric(sub(".*?(\\d+)S$", "\\1", cigar))
  }

  # Return the total length minus the trailing soft clip length
  return(read_len - softclip_len)
}
