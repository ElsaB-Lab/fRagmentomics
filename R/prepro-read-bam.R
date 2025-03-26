#' Preprocess BAM file to extract reads around a given position
#' This function extracts reads from a BAM file using Rsamtools, filtering them
#' based on flags and selecting reads of interest around a genomic position.
#'
#' @param bam A dataframe containing sequencing reads.
#' @param chr Character. Chromosome of interest.
#' @param pos Integer. Genomic position of interest.
#' @param neg_offset Integer. Start position for wider read extraction.
#' @param pos_offset Integer. End position for wider read extraction.
#' @param flag_keep Integer. SAM flag for reads to keep.
#' @param flag_remove Integer. SAM flag for reads to exclude.
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
    neg_offset,
    pos_offset,
    flag_keep,
    flag_remove) {
  # Convert into hexadecimal number
  flag_keep_int <- as.integer(flag_keep)
  flag_remove_int <- as.integer(flag_remove)

  # ---------------------------------------
  # Extract reads around the interested position
  # ---------------------------------------
  region_pos <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = pos, end = pos)
  )
  param_pos <- Rsamtools::ScanBamParam(
    which = region_pos,
    what = c("qname", "flag", "rname", "pos", "mapq", "cigar", "seq", "qual")
  )
  bam_pos_list <- Rsamtools::scanBam(bam, param = param_pos)

  # If we found no reads, we return a empty df
  if (length(bam_pos_list[[1]]$qname) == 0) {
    sam_pos_df <- data.frame(
      QNAME = character(),
      FLAG = integer(),
      RNAME = character(),
      POS = integer(),
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
    (bitwAnd(FLAG, flag_keep_int) == flag_keep_int) &
      (bitwAnd(FLAG, flag_remove_int) == 0)
  )

  # ---------------------------------------
  # Extract wider reads around the interested position
  # ---------------------------------------
  region_ext <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = pos - neg_offset, end = pos + pos_offset)
  )
  param_ext <- Rsamtools::ScanBamParam(
    which = region_ext,
    what = c("qname", "flag", "rname", "pos", "mapq", "cigar", "seq", "qual")
  )
  bam_ext_list <- Rsamtools::scanBam(bam, param = param_ext)

  if (length(bam_ext_list[[1]]$qname) == 0) {
    sam_ext_df <- data.frame(
      QNAME = character(),
      FLAG = integer(),
      RNAME = character(),
      POS = integer(),
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
    (bitwAnd(FLAG, flag_keep_int) == flag_keep_int) &
      (bitwAnd(FLAG, flag_remove_int) == 0)
  )

  # Select only the lines where QNAME (col 1) is in sam_pos_df
  fragments_of_interest <- unique(sam_pos_df$QNAME)
  sam_final_df <- sam_ext_df[sam_ext_df$QNAME %in% fragments_of_interest, ]

  # Return df with sam data
  return(sam_final_df)
}
