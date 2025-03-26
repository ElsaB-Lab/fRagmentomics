# Project : ElsaBLab_fRagmentomics

# Get bit values from bit names of BAM flag.
#'
#' @param flag an integer value
#' @param bitnames an integer vector
#'
#' @return a matrix with column names matching bitnames and binary values
#'
#' @noRd
bitvalues_from_bam_flag <- function(flag, bitnames) {
  FLAG_BITNAMES <- c(
    "isPaired",
    "isProperPair",
    "isUnmappedQuery",
    "hasUnmappedMate",
    "isMinusStrand",
    "isMateMinusStrand",
    "isFirstMateRead",
    "isSecondMateRead",
    "isSecondaryAlignment",
    "isNotPassingQualityControls",
    "isDuplicate",
    "isSupplementaryAlignment"
  )

  bitpos <- match(bitnames, FLAG_BITNAMES)
  invalid_bitnames_idx <- which(is.na(bitpos))
  if (length(invalid_bitnames_idx) != 0L) {
    in1string <- paste0(bitnames[invalid_bitnames_idx], collapse = ", ")
    stop("invalid bitname(s): ", in1string)
  }

  ans <- S4Vectors:::explodeIntBits(flag, bitpos = bitpos)
  dimnames(ans) <- list(names(flag), bitnames)
  ans
}

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

#' Find quality of the base
#'
#' @param current_pos integer representing position before the cigar operations
#' @param r_pos numeric value representing the read mapping position
#' @param cigar character vector representing the read CIGAR
#' @param r_qual character vector representing the read sequencing qualities
#'
#' @return 0 or 1 if the indel is the one we're looking for
#'
#' @noRd
find_c_qual <- function(current_pos, r_pos, r_cigar, r_qual) {
  # Extraction of the sum of the cigar with the S
  matches <- gregexpr("^[0-9]+S", r_cigar) # The ^ means start of the chain
  extracted <- regmatches(r_cigar, matches)[[1]]

  # Put 0 if no soft clipping
  if (length(extracted) == 0) {
    start_soft_clip <- 0
  } else {
    start_soft_clip <- as.numeric(gsub("S", "", extracted))
  }

  pos_indel <- start_soft_clip + current_pos - r_pos

  # Keep the sequence from the nucleotide before the insertion (like in ref)
  insertion_qual <- substr(r_qual, pos_indel, pos_indel)

  insertion_qual
}
