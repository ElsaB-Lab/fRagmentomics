#' Extract and filter paired-end reads for a target locus from a BAM file
#'
#' @description This function performs a targeted extraction of sequencing reads from a BAM file. It first fetches reads
#' within a specified genomic window around an interval of interest, then expands the selection to include the mates of
#' any reads overlapping the interval, ensuring complete fragments are retrieved for analysis.
#'
#' @details
#' The read retrieval and filtering process follows a multi-step pipeline:
#' 1.  A genomic query region is defined around the target interval using the 'neg_offset_mate_search'
#'     and 'pos_offset_mate_search' parameters.
#' 2.  'Rsamtools::scanBam' is used to fetch all reads within this region that pass the filters
#'     specified by 'flag_bam_list'.
#' 3.  From this initial set, the function identifies the subset of reads that *directly* overlap the target interval.
#' 4.  The names ('QNAME') of the fragments corresponding to these overlapping reads are collected.
#' 5.  The function then selects **all** reads from the initially fetched data that share these fragment names, thereby
#'     retrieving the mates even if they did not directly overlap the target interval.
#' 6.  The function returns all reads belonging to the fragments of interest.
#'
#' @inheritParams run_fRagmentomics_variant
#' @param chr Character. Chromosome of interest.
#' @param start Integer. Start of the target interval (1-based, inclusive).
#' @param end Integer. End of the target interval (1-based, inclusive).
#'
#' @return A data.frame of reads belonging to fragments overlapping the target interval,
#'   or \code{NULL} if no reads overlap the interval.
#'
#' @importFrom Rsamtools ScanBamParam scanBam
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @keywords internal
read_bam <- function(
  bam, chr, start, end, neg_offset_mate_search, pos_offset_mate_search,
  flag_bam_list
) {
  # --------------------------------------- Define region and flags for a
  # single, filtered read from the BAM file
  # ---------------------------------------
  start_ext <- max(1, start + neg_offset_mate_search)
  end_ext <- end + pos_offset_mate_search

  region_ext <- GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(
    start = start_ext,
    end = end_ext
  ))

  # Generate the flag vector using the user-provided list
  scan_flag <- do.call(Rsamtools::scanBamFlag, flag_bam_list)

  what_to_scan <- c(
    "qname", "flag", "rname", "pos", "isize", "mapq", "cigar",
    "mrnm", "mpos", "seq", "qual"
  )

  param_ext <- Rsamtools::ScanBamParam(
    which = region_ext, what = what_to_scan,
    flag = scan_flag
  )

  bam_list <- Rsamtools::scanBam(bam, param = param_ext)[[1]]

  if (length(bam_list$qname) == 0) {
    return(NULL)
  }

  # --------------------------------------- Helper function to parse CIGAR
  # strings in base R ---------------------------------------
  get_cigar_width <- function(cigar) {
    matches <- gregexpr("\\d+[MDN=X]", cigar)
    vapply(regmatches(cigar, matches), function(cigar_parts) {
      if (length(cigar_parts) == 0) {
        return(0)
      }
      sum(as.numeric(gsub("[A-Z]", "", cigar_parts)))
    }, FUN.VALUE = numeric(1))
  }

  # --------------------------------------- Convert to a dataframe and find
  # reads covering the position of interest
  # ---------------------------------------
  df_sam_filtered <- data.frame(
    QNAME = as.character(bam_list$qname), FLAG = as.integer(bam_list$flag),
    RNAME = as.character(bam_list$rname), POS = as.integer(bam_list$pos), TLEN = as.integer(bam_list$isize),
    MAPQ = as.integer(bam_list$mapq), CIGAR = as.character(bam_list$cigar), RNEXT = as.character(bam_list$mrnm),
    PNEXT = as.integer(bam_list$mpos), SEQ = as.character(bam_list$seq), QUAL = as.character(bam_list$qual),
    stringsAsFactors = FALSE
  )

  read_width <- get_cigar_width(df_sam_filtered$CIGAR)
  read_end <- df_sam_filtered$POS + read_width - 1

  overlap_interval <- df_sam_filtered$POS <= end & read_end >= start
  fragments_of_interest <- unique(df_sam_filtered[overlap_interval, "QNAME"])

  if (length(fragments_of_interest) == 0) {
    return(NULL)
  }

  # --------------------------------------- Select all reads belonging to the
  # fragments of interest ---------------------------------------
  fragment_name_covering <- df_sam_filtered$QNAME %in% fragments_of_interest
  df_reads_of_covering_fragments <- df_sam_filtered[fragment_name_covering, ]

  return(df_reads_of_covering_fragments)
}
