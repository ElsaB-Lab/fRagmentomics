#' Process a single DNA fragment from paired-end reads in a genomic region
#'
#' @description
#' Worker function that extracts fragment-level features for a fragment
#' overlapping a genomic region. It performs the same quality control,
#' soft-clipping handling, and feature extraction as
#' `extract_fragment_features()`, without any mutation-related information.
#'
#' @inheritParams extract_fragment_features
#' @param region_id Character. Identifier of the region.
#' @param region_start Integer. Start of the region (1-based, inclusive).
#' @param region_end Integer. End of the region (1-based, inclusive).
#'
#' @return A single-row list with the region fragment information.
#'
#' @importFrom Rsamtools bamFlagAsBitMatrix
#'
#' @keywords internal
extract_region_fragment_features <- function(
  df_sam, fragment_name, sample_id, chr, region_id, region_start, region_end,
  report_bam_info, report_softclip, report_5p_3p_bases_fragment, remove_softclip
) {
  mask_frag <- df_sam$QNAME == fragment_name
  df_fragment_reads <- df_sam[mask_frag, , drop = FALSE]

  fragment_qc <- process_fragment_reads_qc(df_fragment_reads, chr)

  if (fragment_qc != "") {
    return(return_fail_qc_region_fragment(
      fragment_qc, sample_id, chr, region_id, region_start, region_end,
      fragment_name, report_bam_info, report_5p_3p_bases_fragment,
      report_softclip
    ))
  }

  flag_matrix <- bamFlagAsBitMatrix(df_fragment_reads$FLAG)
  idx_5p <- which(flag_matrix[, "isMinusStrand"] == 0)
  idx_3p <- which(flag_matrix[, "isMinusStrand"] == 1)

  read_stats_5p <- get_read_stats(df_fragment_reads[idx_5p, ])
  read_stats_3p <- get_read_stats(df_fragment_reads[idx_3p, ])

  if (remove_softclip) {
    input_read_stats_5p <- read_stats_5p
    input_read_stats_3p <- read_stats_3p

    read_5p_info_without_softclip <- remove_softclip(read_stats_5p)
    read_3p_info_without_softclip <- remove_softclip(read_stats_3p)

    if (read_5p_info_without_softclip$SEQ == "" ||
      read_5p_info_without_softclip$CIGAR == "") {
      fragment_qc <- "Invalid read after softclip trimming"
      return(return_fail_qc_region_fragment(
        fragment_qc, sample_id, chr, region_id, region_start, region_end,
        fragment_name, report_bam_info, report_5p_3p_bases_fragment,
        report_softclip
      ))
    }

    read_stats_5p$SEQ <- read_5p_info_without_softclip$SEQ
    read_stats_5p$QUAL <- read_5p_info_without_softclip$QUAL
    read_stats_5p$CIGAR <- read_5p_info_without_softclip$CIGAR
    read_stats_5p$read_length <- read_5p_info_without_softclip$read_length

    if (read_3p_info_without_softclip$SEQ == "" ||
      read_3p_info_without_softclip$CIGAR == "") {
      fragment_qc <- "Invalid read after softclip trimming"
      return(return_fail_qc_region_fragment(
        fragment_qc, sample_id, chr, region_id, region_start, region_end,
        fragment_name, report_bam_info, report_5p_3p_bases_fragment,
        report_softclip
      ))
    }

    read_stats_3p$SEQ <- read_3p_info_without_softclip$SEQ
    read_stats_3p$QUAL <- read_3p_info_without_softclip$QUAL
    read_stats_3p$CIGAR <- read_3p_info_without_softclip$CIGAR
    read_stats_3p$read_length <- read_3p_info_without_softclip$read_length
  }

  fragment_size <- get_fragment_size(read_stats_5p, read_stats_3p)
  Position_3p <- end_on_reference(read_stats_3p$POS, read_stats_3p$CIGAR)

  output_read_stats_5p <- if (remove_softclip) input_read_stats_5p else read_stats_5p
  output_read_stats_3p <- if (remove_softclip) input_read_stats_3p else read_stats_3p

  final_row_fragment <- return_fail_qc_region_fragment(
    "OK", sample_id, chr, region_id, region_start, region_end,
    fragment_name, report_bam_info, report_5p_3p_bases_fragment,
    report_softclip
  )
  final_row_fragment$Fragment_Size <- as.integer(fragment_size)
  final_row_fragment$Position_5p <- as.integer(output_read_stats_5p$POS)
  final_row_fragment$Position_3p <- as.integer(Position_3p)

  if (report_bam_info) {
    final_row_fragment$POS_5p <- as.integer(output_read_stats_5p$POS)
    final_row_fragment$POS_3p <- as.integer(output_read_stats_3p$POS)
    final_row_fragment$FLAG_5p <- as.integer(output_read_stats_5p$FLAG)
    final_row_fragment$FLAG_3p <- as.integer(output_read_stats_3p$FLAG)
    final_row_fragment$MAPQ_5p <- as.integer(output_read_stats_5p$MAPQ)
    final_row_fragment$MAPQ_3p <- as.integer(output_read_stats_3p$MAPQ)
    final_row_fragment$CIGAR_5p <- as.character(output_read_stats_5p$CIGAR)
    final_row_fragment$CIGAR_3p <- as.character(output_read_stats_3p$CIGAR)

    if (is.null(read_stats_5p$TLEN) || is.na(read_stats_5p$TLEN)) {
      warning("TLEN is NULL or missing, cannot be reported.")
      final_row_fragment$TLEN <- NA_integer_
    } else {
      final_row_fragment$TLEN <- as.integer(abs(read_stats_5p$TLEN))
    }
  }

  if (report_5p_3p_bases_fragment != 0) {
    fragment_bases_5p_3p <- get_fragment_bases_5p_3p(
      report_5p_3p_bases_fragment,
      read_stats_5p$SEQ,
      read_stats_3p$SEQ,
      read_stats_5p$QUAL,
      read_stats_3p$QUAL
    )

    final_row_fragment$Fragment_Bases_5p <-
      as.character(fragment_bases_5p_3p$fragment_bases_5p)
    final_row_fragment$Fragment_Bases_3p <-
      as.character(fragment_bases_5p_3p$fragment_bases_3p)
    final_row_fragment$Fragment_Basqs_5p <-
      as.character(fragment_bases_5p_3p$fragment_basqs_5p)
    final_row_fragment$Fragment_Basqs_3p <-
      as.character(fragment_bases_5p_3p$fragment_basqs_3p)
  }

  if (report_softclip) {
    fragment_bases_softclip_5p_3p <- get_fragment_bases_5p_3p_softclip(
      read_stats_5p$CIGAR,
      read_stats_3p$CIGAR
    )

    final_row_fragment$Nb_Fragment_Bases_Softclip_5p <-
      as.integer(fragment_bases_softclip_5p_3p$nb_softclip_5p)
    final_row_fragment$Nb_Fragment_Bases_Softclip_3p <-
      as.integer(fragment_bases_softclip_5p_3p$nb_softclip_3p)
  }

  return(final_row_fragment)
}

#' Generate a standardized failure result for a region fragment analysis
#'
#' @description
#' Creates the region fragment output schema with `NA` values for a fragment
#' that failed quality control. Optional columns follow the requested reporting
#' flags so that failed and successful rows share the same schema.
#'
#' @param qc_message Character. Quality control message (`"OK"` for success).
#' @inheritParams extract_region_fragment_features
#'
#' @return A named list with the region fragment schema.
#'
#' @noRd
return_fail_qc_region_fragment <- function(
  qc_message, sample_id, chr, region_id, region_start, region_end,
  fragment_name, report_bam_info, report_5p_3p_bases_fragment, report_softclip
) {
  row <- list(
    Sample_Id = if (is.na(sample_id)) NA_character_ else as.character(sample_id),
    Chromosome = chr,
    Region_Id = region_id,
    Region_Start = as.integer(region_start),
    Region_End = as.integer(region_end),
    Region_Width = as.integer(region_end - region_start + 1L),
    Fragment_Id = fragment_name,
    Fragment_QC = qc_message,
    Fragment_Size = NA_integer_,
    Position_5p = NA_integer_,
    Position_3p = NA_integer_
  )

  if (report_bam_info) {
    row$POS_5p <- NA_integer_
    row$POS_3p <- NA_integer_
    row$FLAG_5p <- NA_integer_
    row$FLAG_3p <- NA_integer_
    row$MAPQ_5p <- NA_integer_
    row$MAPQ_3p <- NA_integer_
    row$CIGAR_5p <- NA_character_
    row$CIGAR_3p <- NA_character_
    row$TLEN <- NA_integer_
  }

  if (report_5p_3p_bases_fragment != 0) {
    row$Fragment_Bases_5p <- NA_character_
    row$Fragment_Bases_3p <- NA_character_
    row$Fragment_Basqs_5p <- NA_character_
    row$Fragment_Basqs_3p <- NA_character_
  }

  if (report_softclip) {
    row$Nb_Fragment_Bases_Softclip_5p <- NA_integer_
    row$Nb_Fragment_Bases_Softclip_3p <- NA_integer_
  }

  row
}
