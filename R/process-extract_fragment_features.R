#' Process a single sequencing fragment
#' Extracts and processes relevant data from a sequencing fragment.
#'
#' @inheritParams analyze_fragments
#' @param df_sam A dataframe containing sequencing reads.
#' @param fragment_name Name of the fragment (paired-end reads).
#' @param chr Character vector representing the chromosome of interest.
#' @param pos Numeric value representing the Genomic position of interest.
#' @param ref Character vector representing reference base(s).
#' @param alt Character vector representing alternative base(s).
#' @param fasta_fafile An open connection to an object of class FaFile
#'
#' @return A dataframe with the processed fragment information.
#'
#' @export
extract_fragment_features <- function(df_sam,
                                      fragment_name,
                                      sample_id,
                                      chr,
                                      pos,
                                      ref,
                                      alt,
                                      report_tlen,
                                      report_softclip,
                                      report_5p_3p_bases_fragment,
                                      cigar_free_mode,
                                      fasta_fafile) {
  # Select reads originating from the fragment of interest
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

  # If the fragment fails QC, return a dataframe with NAs
  if (fragment_qc != "") {
    final_row_fragment <- list(
      Chromosome = chr,
      Position = pos,
      Ref = ref,
      Alt = alt,
      Fragment_Id = fragment_name,
      Fragment_QC = fragment_qc,
      Fragment_Status = NA,
      Fragment_Status_Broad = NA,
      Fragment_size = NA,
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

  # Separate read_1/read_2
  reads <- identify_read_sequence_order(df_fragment_reads)
  read_stats_1 <- get_read_stats(reads$read_1)
  read_stats_2 <- get_read_stats(reads$read_2)

  # -------------------------------
  # Get read sequence, read base qualities, and read mutation status
  # -------------------------------
  read_info_1 <- get_base_basq_mstat_from_read(chr, pos, ref, alt, read_stats_1, fasta_fafile, cigar_free_mode)
  read_info_2 <- get_base_basq_mstat_from_read(chr, pos, ref, alt, read_stats_2, fasta_fafile, cigar_free_mode)

  # -------------------------------
  # Compute fragment size
  # -------------------------------
  inner_distance <- get_fragment_inner_distance(
    read_stats_1$POS, read_stats_1$CIGAR, read_stats_1$read_length,
    read_stats_2$POS, read_stats_2$CIGAR, read_stats_2$read_length
  )
  absolute_size <- read_stats_1$read_length +
    read_stats_2$read_length +
    inner_distance

  # -------------------------------
  # Define fragment status
  # -------------------------------
  fragment_status_broad <- get_fragment_mutation_status_broad(mstat_1 = read_info_1$mstat, mstat_2 = read_info_2$mstat)
  fragment_status <- get_fragment_mutation_status(fragment_status_broad)

  # -------------------------------
  # Define 3' and 5' reads
  # -------------------------------
  identity_3p <- ifelse(read_stats_1$POS <= read_stats_2$POS, 2, 1)
  identity_5p <- ifelse(identity_3p == 2, 1, 2)

  if (identity_5p == 1) {
    read_stats_5p <- read_stats_1
    read_info_5p <- read_info_1
    read_stats_3p <- read_stats_2
    read_info_3p <- read_info_2
  } else {
    read_stats_5p <- read_stats_2
    read_info_5p <- read_info_2
    read_stats_3p <- read_stats_1
    read_info_3p <- read_info_1
  }

  # -------------------------------
  # Build an adaptative dataframe
  # -------------------------------
  final_row_fragment <- list(
    Chromosome            = chr,
    Position              = pos,
    Ref                   = ref,
    Alt                   = alt,
    Fragment_Id           = fragment_name,
    Fragment_QC           = "OK",
    Fragment_Status       = fragment_status,
    Fragment_Status_Broad = fragment_status_broad,
    Fragment_Size         = absolute_size,
    Inner_Distance        = inner_distance,
    Read_5p               = identity_5p,
    Read_5p_Status        = read_info_5p$mstat,
    Read_3p_Status        = read_info_3p$mstat,
    MAPQ_5p               = read_stats_5p$MAPQ,
    MAPQ_3p               = read_stats_3p$MAPQ,
    BASE_5p               = read_info_5p$base,
    BASE_3p               = read_info_3p$base,
    BASQ_5p               = read_info_5p$basq,
    BASQ_3p               = read_info_3p$basq,
    CIGAR_5p              = read_stats_5p$CIGAR,
    CIGAR_3p              = read_stats_3p$CIGAR,
    POS_5p                = read_stats_5p$POS,
    POS_3p                = read_stats_3p$POS
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
    if (is.null(read_stats_1$TLEN) || is.na(read_stats_1$TLEN)) {
      message("Warning: TLEN is NULL")
      final_row_fragment$TLEN <- "Warning: TLEN is NULL"
    } else {
      final_row_fragment$TLEN <- abs(read_stats_1$TLEN)
    }
  }

  # -------------------------------
  # Define n base 5' and 3'
  # -------------------------------
  # if report_5p_3p_bases_fragment != 0,
  # Call the function get_fragment_bases_5p_3p
  if (report_5p_3p_bases_fragment != 0) {
    fragment_bases_5p_3p <- get_fragment_bases_5p_3p(
      report_5p_3p_bases_fragment,
      read_stats_5p$SEQ,
      read_stats_3p$SEQ,
      read_stats_5p$QUAL,
      read_stats_3p$QUAL
    )

    final_row_fragment$Fragment_Bases_5p <-
      fragment_bases_5p_3p$fragment_bases_5p
    final_row_fragment$Fragment_Bases_3p <-
      fragment_bases_5p_3p$fragment_bases_3p
    final_row_fragment$Fragment_Basqs_5p <-
      fragment_bases_5p_3p$fragment_basqs_5p
    final_row_fragment$Fragment_Basqs_3p <-
      fragment_bases_5p_3p$fragment_basqs_3p
  }

  # -------------------------------
  # Define number of soft clipped bases in 5' and 3'
  # -------------------------------
  # If report_softclip is TRUE,
  # Call the function get_fragment_bases_5p_3p_softclip
  if (report_softclip) {
    fragment_bases_softclip_5p_3p <- get_fragment_bases_5p_3p_softclip(
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


#' Identify the first read and second read order using the BAM FLAG field.
#'
#' @param df_fragment_reads A dataframe with two rows representing reads.
#'
#' @return A named list containing two elements: `read_1` and `read_2`.
#'
#' @noRd
identify_read_sequence_order <- function(df_fragment_reads) {
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
    list(read_1 = df_fragment_reads[1, ], read_2 = df_fragment_reads[2, ])
  } else if (
    (flag_b[, "isFirstMateRead"] == 1) &&
      (flag_a[, "isSecondMateRead"] == 1)
  ) {
    list(read_1 = df_fragment_reads[2, ], read_2 = df_fragment_reads[1, ])
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
#' mapq, tlen, cigar, pos, seq, qual, and read_length.
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
