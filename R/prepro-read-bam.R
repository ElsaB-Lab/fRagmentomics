# Project : ElsaBLab_fRagmentomics

#' Preprocess BAM file to extract reads around a given position
#' This function extracts reads from a BAM file using Rsamtools, filtering them
#' based on flags and selecting reads of interest around a genomic position.
#'
#' @inheritParams process_fragmentomics
#' @param chr Character. Chromosome of interest.
#' @param pos Integer. Genomic position of interest.
#'
#' @return A data frame containing the filtered SAM entries.
#'
#' @importFrom Rsamtools ScanBamParam
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
#' @keywords internal
read_bam <- function(
    bam,
    chr,
    pos,
    neg_offset_mate_search,
    pos_offset_mate_search,
    flag_keep,
    flag_remove) {
  # Convert into hexadecimal number
  flag_keep_int <- as.integer(flag_keep)
  flag_remove_int <- as.integer(flag_remove)

  # ---------------------------------------
  # Check the start and end position of the read bam
  # ---------------------------------------
  start <- max(1, pos + neg_offset_mate_search)

  # ---------------------------------------
  # Extract reads around the interested position
  # ---------------------------------------
  region_pos <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = pos, end = pos)
  )
  param_pos <- Rsamtools::ScanBamParam(
    which = region_pos,
    what = c("qname", "flag", "rname", "pos", "isize", "mapq", "cigar", "seq", "qual")
  )
  bam_pos_list <- Rsamtools::scanBam(bam, param = param_pos)

  # If we found no reads, we return a empty df
  if (length(bam_pos_list[[1]]$qname) == 0) {
    sam_pos_df <- data.frame(
      QNAME = character(),
      FLAG = integer(),
      RNAME = character(),
      POS = integer(),
      TLEN = integer(),
      MAPQ = integer(),
      CIGAR = character(),
      SEQ = character(),
      QUAL = character(),
      stringsAsFactors = FALSE
    )
  } else {
    sam_pos_df <- data.frame(
      QNAME = as.character(bam_pos_list[[1]]$qname),
      FLAG = as.integer(bam_pos_list[[1]]$flag),
      RNAME = as.character(bam_pos_list[[1]]$rname),
      POS = as.integer(bam_pos_list[[1]]$pos),
      TLEN = as.integer(bam_pos_list[[1]]$isize),
      MAPQ = as.integer(bam_pos_list[[1]]$mapq),
      CIGAR = as.character(bam_pos_list[[1]]$cigar),
      SEQ = as.character(bam_pos_list[[1]]$seq),
      QUAL = as.character(bam_pos_list[[1]]$qual),
      stringsAsFactors = FALSE
    )
  }

  # We only keep the reads with the appropriate flags
  sam_pos_df <- subset(
    sam_pos_df,
    (bitwAnd(sam_pos_df$FLAG, flag_keep_int) == flag_keep_int) &
      (bitwAnd(sam_pos_df$FLAG, flag_remove_int) == 0)
  )

  # ---------------------------------------
  # Extract wider reads around the interested position
  # ---------------------------------------
  region_ext <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = start, end = pos + pos_offset_mate_search)
  )
  param_ext <- Rsamtools::ScanBamParam(
    which = region_ext,
    what = c("qname", "flag", "rname", "pos", "isize", "mapq", "cigar", "seq", "qual")
  )
  bam_ext_list <- Rsamtools::scanBam(bam, param = param_ext)

  if (length(bam_ext_list[[1]]$qname) == 0) {
    sam_ext_df <- data.frame(
      QNAME = character(),
      FLAG = integer(),
      RNAME = character(),
      POS = integer(),
      TLEN = integer(),
      MAPQ = integer(),
      CIGAR = character(),
      SEQ = character(),
      QUAL = character(),
      stringsAsFactors = FALSE
    )
  } else {
    sam_ext_df <- data.frame(
      QNAME = as.character(bam_ext_list[[1]]$qname),
      FLAG = as.integer(bam_ext_list[[1]]$flag),
      RNAME = as.character(bam_ext_list[[1]]$rname),
      POS = as.integer(bam_ext_list[[1]]$pos),
      TLEN = as.integer(bam_ext_list[[1]]$isize),
      MAPQ = as.integer(bam_ext_list[[1]]$mapq),
      CIGAR = as.character(bam_ext_list[[1]]$cigar),
      SEQ = as.character(bam_ext_list[[1]]$seq),
      QUAL = as.character(bam_ext_list[[1]]$qual),
      stringsAsFactors = FALSE
    )
  }

  # We only keep the reads with the appropriate flags
  sam_ext_df <- subset(
    sam_ext_df,
    (bitwAnd(sam_ext_df$FLAG, flag_keep_int) == flag_keep_int) &
      (bitwAnd(sam_ext_df$FLAG, flag_remove_int) == 0)
  )

  # Select only the lines where QNAME (col 1) is in sam_pos_df
  fragments_of_interest <- unique(sam_pos_df$QNAME)
  sam_final_df <- sam_ext_df[sam_ext_df$QNAME %in% fragments_of_interest, ]

  # Check if the final bam is empty and have the selected columns
  if (nrow(sam_final_df) == 0) {
    stop("The final BAM data frame is empty. No reads match the criteria.")
  }

  required_columns <-
    c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "SEQ", "QUAL")
  missing_columns <- setdiff(required_columns, colnames(sam_final_df))

  if (length(missing_columns) > 0) {
    stop(paste(
      "The final BAM data frame is missing the following columns:",
      paste(missing_columns, collapse = ", ")
    ))
  }

  # Return df with sam data
  sam_final_df
}
