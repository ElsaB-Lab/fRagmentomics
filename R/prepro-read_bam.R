#' Preprocess BAM file to extract reads around a given position
#' This function extracts reads from a BAM file using Rsamtools, filtering them
#' based on flags and selecting reads of interest around a genomic position.
#'
#' @inheritParams analyze_fragments
#' @param chr Character. Chromosome of interest.
#' @param pos Integer. Genomic position of interest.
#'
#' @return A dataframe containing the filtered SAM entries.
#'
#' @importFrom Rsamtools ScanBamParam scanBam
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
    flag_bam_list) {
  # ---------------------------------------
  # Define region and flags for a single, filtered read from the BAM file
  # ---------------------------------------
  start_ext <- max(1, pos + neg_offset_mate_search)
  end_ext <- pos + pos_offset_mate_search

  region_ext <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = start_ext, end = end_ext)
  )

  # Generate the flag vector using the user-provided list
  scan_flag <- do.call(Rsamtools::scanBamFlag, flag_bam_list)

  what_to_scan <- c("qname", "flag", "rname", "pos", "isize", "mapq", "cigar", "seq", "qual")

  param_ext <- Rsamtools::ScanBamParam(
    which = region_ext,
    what = what_to_scan,
    flag = scan_flag
  )

  bam_list <- Rsamtools::scanBam(bam, param = param_ext)[[1]]

  if (length(bam_list$qname) == 0) {
      return(NULL)
  }

  # ---------------------------------------
  # Helper function to parse CIGAR strings in base R
  # ---------------------------------------
  get_cigar_width <- function(cigar) {
    matches <- gregexpr("\\d+[MDN=X]", cigar)
    sapply(regmatches(cigar, matches), function(cigar_parts) {
      if (length(cigar_parts) == 0) {
        return(0)
      }
      sum(as.numeric(gsub("[A-Z]", "", cigar_parts)))
    })
  }

  # ---------------------------------------
  # Convert to a dataframe and find reads covering the position of interest
  # ---------------------------------------
  df_sam_filtered <- data.frame(
    QNAME = as.character(bam_list$qname),
    FLAG = as.integer(bam_list$flag),
    RNAME = as.character(bam_list$rname), 
    POS = as.integer(bam_list$pos),
    TLEN = as.integer(bam_list$isize),
    MAPQ = as.integer(bam_list$mapq),
    CIGAR = as.character(bam_list$cigar),
    SEQ = as.character(bam_list$seq),
    QUAL = as.character(bam_list$qual),
    stringsAsFactors = FALSE
  )

  read_width <- get_cigar_width(df_sam_filtered$CIGAR)
  read_end <- df_sam_filtered$POS + read_width - 1

  cover_position <- df_sam_filtered$POS <= pos & read_end >= pos
  fragments_of_interest <- unique(df_sam_filtered[cover_position, "QNAME"])

  if (length(fragments_of_interest) == 0) {
    return(NULL) 
  }

  # ---------------------------------------
  # Select all reads belonging to the fragments of interest
  # ---------------------------------------
  final_condition <- df_sam_filtered$QNAME %in% fragments_of_interest
  df_sam_final <- df_sam_filtered[final_condition, ]

  return(df_sam_final)
}
