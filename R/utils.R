# Project : ElsaBLab_fRagmentomics

#' @title bitvalues_from_bam_flag
#' @description Get bit values from bit names of BAM flag.
#'
#' @param flag an integer value
#' @param bitnames an integer vector
#' @return a matrix with column names matching bitnames and binary values representing bit values
#'
#' @keywords internal
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
      in1string <- paste0(bitnames[invalid_bitnames_idx], collapse=", ")
      stop("invalid bitname(s): ", in1string)
  }

  ans <- S4Vectors:::explodeIntBits(flag, bitpos=bitpos)
  dimnames(ans) <- list(names(flag), bitnames)
  ans
}

#' @title parse_cigar
#' @description Function to parse the CIGAR string
#'
#' @param cigar a character vector
#' @return a df with length and type of the CIGAR string. One row per operation.
#'
#' @keywords internal
parse_cigar <- function(cigar) {
  # Extract the operations and lengths from the CIGAR string
  matches <- gregexpr("([0-9]+)([MIDNSHP=X])", cigar, perl=TRUE)
  ops_str <- regmatches(cigar, matches)[[1]]

  lengths <- as.numeric(gsub("([0-9]+)([MIDNSHP=X])", "\\1", ops_str))
  types <- gsub("([0-9]+)([MIDNSHP=X])", "\\2", ops_str)

  data.frame(length = lengths, type = types, stringsAsFactors = FALSE)
}


