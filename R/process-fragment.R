# Project : ElsaBLab_fRagmentomics

#' Identify read1 and read2 from flags
#' Determines which read is the first and seconde mate (read1 and read2).
#'
#' @param df_fragment_reads A dataframe with two rows representing reads.
#'
#' @return A named list containing two elements: `read1` and `read2`.
#'
#' @noRd
identify_reads <- function(df_fragment_reads) {
  # Extract flag values
  flag_a <- bitvalues_from_bam_flag(
    as.integer(df_fragment_reads[1, "FLAG"]),
    bitnames = c("isFirstMateRead", "isSecondMateRead")
  )
  flag_b <- bitvalues_from_bam_flag(
    as.integer(df_fragment_reads[2, "FLAG"]),
    bitnames = c("isFirstMateRead", "isSecondMateRead")
  )

  if (
    (flag_a[, "isFirstMateRead"] == 1) &&
      (flag_b[, "isSecondMateRead"] == 1)
  ) {
    list(read1 = df_fragment_reads[1, ], read2 = df_fragment_reads[2, ])
  } else if (
    (flag_b[, "isFirstMateRead"] == 1) &&
      (flag_a[, "isSecondMateRead"] == 1)
  ) {
    list(read1 = df_fragment_reads[2, ], read2 = df_fragment_reads[1, ])
  } else {
    stop("Invalid read flags.
    One read should be first mate, the other second mate.")
  }
}

#' Retrieves key attributes from a single sequencing read.
#'
#' @param df_read A dataframe containing a single read (one row).
#'
#' @return A list with extracted read properties:
#' mapq, tlen, cigar, pos, query, qual, and read_length.
#'
#' @noRd
get_read_stats <- function(df_read) {
  list(
    MAPQ        = df_read$MAPQ,
    TLEN        = df_read$TLEN,
    CIGAR       = df_read$CIGAR,
    POS         = df_read$POS,
    SEQ         = df_read$SEQ,
    QUAL        = df_read$QUAL,
    read_length = nchar(df_read$SEQ)
  )
}

#' Process a single sequencing fragment
#' Extracts and processes relevant data from a sequencing fragment.
#'
#' @inheritParams process_fragmentomics
#' @param df_sam A dataframe containing sequencing reads.
#' @param fragment_name Name of the fragment (paired-end reads).
#' @param chr Character vector representing the chromosome of interest.
#' @param pos Numeric value representing the Genomic position of interest.
#' @param ref Character vector representing reference base(s).
#' @param alt Character vector representing alternative base(s).
#' @param mutation_type In "SNV", "deletion", or "insertion".
#'
#' @return A dataframe with the processed fragment information.
#'
#' @export
process_fragment <- function(df_sam,
                             fragment_name,
                             sample_id,
                             chr,
                             pos,
                             ref,
                             alt,
                             mutation_type,
                             report_tlen,
                             report_softclip,
                             report_5p_3p_bases_fragment) {
  df_fragment_reads <- df_sam[
    df_sam[, 1, drop = TRUE] == fragment_name, ,
    drop = FALSE
  ]

  # -------------------------------
  # Sanity check fragments
  # -------------------------------
  fragment_qc <- ""
  if (!all(df_fragment_reads[, 3] == chr)) {
    fragment_qc <- paste(fragment_qc, "& Read(s) from another chromosome
    than", chr)
  }
  if (nrow(df_fragment_reads) != 2) {
    fragment_qc <- paste(fragment_qc, "&", nrow(df_fragment_reads), "read(s)")
  }

  if (fragment_qc != "") {
    final_row_fragment <- list(
      Chromosome = chr,
      Position = pos,
      Ref = ref,
      Alt = alt,
      Fragment_Id = fragment_name,
      Fragment_QC = fragment_qc,
      Fragment_Status = NA,
      Absolute_size = NA,
      Inner_distance = NA,
      Read_5p = NA,
      MAPQ_5p = NA,
      MAPQ_3p = NA,
      BASE_5p = NA,
      BASE_3p = NA,
      BASQ_5p = NA,
      BASQ_3p = NA,
      CIGAR_5p = NA,
      CIGAR_3p = NA,
      Pos_bam_5p = NA,
      Pos_bam_3p = NA
    )

    if (!is.na(sample_id)) {
      final_row_fragment$Sample_Id <- sample_id
    }

    if (report_tlen) {
      final_row_fragment$TLEN <- NA
    }

    if (report_5p_3p_bases_fragment != 0) {
      final_row_fragment$Fragment_bases_5p <- NA
      final_row_fragment$Fragment_bases_3p <- NA
      final_row_fragment$Fragment_Qbases_5p <- NA
      final_row_fragment$Fragment_Qbases_3p <- NA
    }

    if (report_softclip) {
      final_row_fragment$Nb_fragment_bases_softclip_5p <- NA
      final_row_fragment$Nb_fragment_bases_softclip_3p <- NA
    }

    result_df <- as.data.frame(final_row_fragment)
    return(result_df)
  }

  # Separate read1/read2
  reads <- identify_reads(df_fragment_reads)
  read1_stats <- get_read_stats(reads$read1)
  read2_stats <- get_read_stats(reads$read2)

  # -------------------------------
  # Retrieve mutation information (a list with base and qual)
  # -------------------------------
  r_info1 <- get_mutation_info(mutation_type, pos, ref, alt, read1_stats)
  r_info2 <- get_mutation_info(mutation_type, pos, ref, alt, read2_stats)

  # -------------------------------
  # Compute fragment size
  # -------------------------------
  inner_distance <- process_fragment_inner_distance(
    read1_stats$POS, read1_stats$CIGAR, read1_stats$read_length,
    read2_stats$POS, read2_stats$CIGAR, read2_stats$read_length
  )
  absolute_size <- read1_stats$read_length +
    read2_stats$read_length +
    inner_distance

  # -------------------------------
  # Define fragment status
  # -------------------------------
  fragment_status <- process_fragment_status(
    ref,
    alt,
    mutation_type,
    r_info1$base,
    r_info2$base
  )

  # -------------------------------
  # Define 3' and 5' reads
  # -------------------------------
  identity_3p <- ifelse(read1_stats$POS <= read2_stats$POS, 2, 1)
  identity_5p <- ifelse(identity_3p == 2, 1, 2)

  if (identity_5p == 1) {
    read_stats_5p <- read1_stats
    r_info_5p <- r_info1
    read_stats_3p <- read2_stats
    r_info_3p <- r_info2
  } else {
    read_stats_5p <- read2_stats
    r_info_5p <- r_info2
    read_stats_3p <- read1_stats
    r_info_3p <- r_info1
  }

  # -------------------------------
  # Build an adaptative dataframe
  # -------------------------------
  final_row_fragment <- list(
    Chromosome       = chr,
    Position         = pos,
    Ref              = ref,
    Alt              = alt,
    Fragment_Id      = fragment_name,
    Fragment_QC      = "OK",
    Fragment_Status  = fragment_status,
    Absolute_size    = absolute_size,
    Inner_distance   = inner_distance,
    Read_5p          = identity_5p,
    MAPQ_5p          = read_stats_5p$MAPQ,
    MAPQ_3p          = read_stats_3p$MAPQ,
    BASE_5p          = r_info_5p$base,
    BASE_3p          = r_info_3p$base,
    BASQ_5p          = r_info_5p$qual,
    BASQ_3p          = r_info_3p$qual,
    CIGAR_5p         = read_stats_5p$CIGAR,
    CIGAR_3p         = read_stats_3p$CIGAR,
    Pos_bam_5p       = read_stats_5p$POS,
    Pos_bam_3p       = read_stats_3p$POS
  )

  # -------------------------------
  # Put sample if not NA
  # -------------------------------
  if (!is.na(sample_id)) {
    final_row_fragment$Sample_Id <- sample_id
  }

  # -------------------------------
  # Put sample if not NA
  # -------------------------------
  if (report_tlen) {
    if (is.null(read1_stats$TLEN) || is.na(read1_stats$TLEN)) {
      message("Warning: TLEN is NULL")
      final_row_fragment$TLEN <- "Warning: TLEN is NULL"
    } else {
      final_row_fragment$TLEN <- abs(read1_stats$TLEN)
    }
  }

  # -------------------------------
  # Define n base 5' and 3'
  # -------------------------------
  # if report_5p_3p_bases_fragment != 0,
  # Call the function process_get_fragment_bases
  if (report_5p_3p_bases_fragment != 0) {
    fragment_bases_5p_3p <- process_get_fragment_bases(
      report_5p_3p_bases_fragment,
      read_stats_5p$SEQ,
      read_stats_3p$SEQ,
      read_stats_5p$QUAL,
      read_stats_3p$QUAL
    )

    final_row_fragment$Fragment_bases_5p <-
      fragment_bases_5p_3p$fragment_bases_5p
    final_row_fragment$Fragment_bases_3p <-
      fragment_bases_5p_3p$fragment_bases_3p
    final_row_fragment$Fragment_Qbases_5p <-
      fragment_bases_5p_3p$fragment_Qbases_5p
    final_row_fragment$Fragment_Qbases_3p <-
      fragment_bases_5p_3p$fragment_Qbases_3p
  }

  # -------------------------------
  # Define number of soft clipped bases in 5' and 3'
  # -------------------------------
  # If report_softclip is TRUE,
  # Call the function process_get_fragment_bases_softclip
  if (report_softclip) {
    fragment_bases_softclip_5p_3p <- process_get_fragment_bases_softclip(
      read_stats_5p$CIGAR,
      read_stats_3p$CIGAR
    )

    final_row_fragment$Nb_fragment_bases_softclip_5p <-
      fragment_bases_softclip_5p_3p$nb_softclip_5p
    final_row_fragment$Nb_fragment_bases_softclip_3p <-
      fragment_bases_softclip_5p_3p$nb_softclip_3p
  }

  # Creation of the dataframe with the wanted columns
  result_df <- as.data.frame(final_row_fragment)
  result_df
}
