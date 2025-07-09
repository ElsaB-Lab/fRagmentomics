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
#' @param fasta_seq A list with the fasta sequence between two positions.
#' @param input_mutation_info Character vecotr representing the input mutation.
#'
#' @return A dataframe with the processed fragment information.
#'
#' @keywords internal
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
                                      cigar_free_indel_match,
                                      fasta_fafile = NULL,
                                      fasta_seq = NULL,
                                      input_mutation_info) {
  # Select reads originating from the fragment of interest
  df_fragment_reads <- df_sam[
    df_sam[, 1, drop = TRUE] == fragment_name, ,
    drop = FALSE
  ]

  # -------------------------------
  # Sanity check fragments
  # -------------------------------
  qc_messages <- list()

  # Test 1: Report if the fragment does not consist of exactly two reads.
  if (nrow(df_fragment_reads) != 2) {
    qc_messages <- c(qc_messages, paste("Fragment has", nrow(df_fragment_reads), "read(s)"))
  }

  # Test 2: Report if any read in the fragment is not on the expected chromosome.
  if (!all(df_fragment_reads$RNAME == chr)) {
    qc_messages <- c(qc_messages, paste("Read(s) found on a chromosome other than", chr))
  }

  if (nrow(df_fragment_reads) > 0) {
    # Isolate the first read of the fragment
    read <- df_fragment_reads[1, , drop = FALSE]

    # Test 3: Report if the mate appears to be on another chromosome (translocation).
    if (!is.na(read$RNEXT) && read$RNEXT != "=" && read$RNAME != read$RNEXT) {
      qc_messages <- c(qc_messages, paste("Mate maps to a different chromosome. RNEXT=", read$RNEXT))
    }

    # Test 4: Report if the read or its mate is marked as unmapped.
    if (read$POS == 0 || read$PNEXT == 0) {
      qc_messages <- c(qc_messages, paste("Read or mate is unmapped. POS=", read$POS, "& PNEXT=", read$PNEXT))
    }
  }

  fragment_qc <- paste(qc_messages, collapse = " & ")

  # If the fragment fails QC, return a dataframe with NAs
  if (fragment_qc != "") {
    final_row_fragment <- list(
      Chromosome             = chr,
      Position               = pos,
      Ref                    = ref,
      Alt                    = alt,
      Input_Mutation         = input_mutation_info,
      Fragment_Id            = fragment_name,
      Fragment_QC            = fragment_qc,
      Fragment_Status_Simple = NA,
      Fragment_Status_Detail = NA,
      Fragment_Size          = NA,
      Inner_Distance         = NA,
      Read_5p                = NA,
      Read_5p_Status         = NA,
      Read_3p_Status         = NA,
      MAPQ_5p                = NA,
      MAPQ_3p                = NA,
      BASE_5p                = NA,
      BASE_3p                = NA,
      BASQ_5p                = NA,
      BASQ_3p                = NA,
      CIGAR_5p               = NA,
      CIGAR_3p               = NA,
      POS_5p                 = NA,
      POS_3p                 = NA
    )

    if (!is.na(sample_id)) {
      final_row_fragment$Sample_Id <- sample_id
    }

    if (report_tlen) {
      final_row_fragment$TLEN <- NA
    }

    if (report_5p_3p_bases_fragment != 0) {
      final_row_fragment$Fragment_Bases_5p <- NA
      final_row_fragment$Fragment_Bases_3p <- NA
      final_row_fragment$Fragment_Basqs_5p <- NA
      final_row_fragment$Fragment_Basqs_3p <- NA
    }

    if (report_softclip) {
      final_row_fragment$Nb_Fragment_Bases_Softclip_5p <- NA
      final_row_fragment$Nb_Fragment_Bases_Softclip_3p <- NA
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
  read_info_1 <- get_base_basq_mstat_from_read(
    chr, pos, ref, alt, read_stats_1, fasta_fafile, fasta_seq,
    cigar_free_indel_match
  )
  read_info_2 <- get_base_basq_mstat_from_read(
    chr, pos, ref, alt, read_stats_2, fasta_fafile, fasta_seq,
    cigar_free_indel_match
  )

  # -------------------------------
  # Define 3' and 5' reads
  # -------------------------------
  # 5p read: forward strand / 3p: reverse strand

  # If bitwAnd(FLAG, 16) == 0, strand is '+'.
  is_read1_on_forward_strand <- bitwAnd(read_stats_1$FLAG, 16) == 0

  if (is_read1_on_forward_strand) {
    # read_1 is on strand '+'. read_1 -> 5p
    read_stats_5p <- read_stats_1
    read_info_5p <- read_info_1
    read_stats_3p <- read_stats_2
    read_info_3p <- read_info_2
  } else {
    # read_1 is on strand '-'. read_1 -> 3p
    read_stats_5p <- read_stats_2
    read_info_5p <- read_info_2
    read_stats_3p <- read_stats_1
    read_info_3p <- read_info_1
  }

  # -------------------------------
  # Compute fragment size
  # -------------------------------
  inner_distance <- get_fragment_inner_distance(
    read_stats_5p$POS, read_stats_5p$CIGAR,
    read_stats_3p$POS, read_stats_3p$CIGAR
  )
  absolute_size <- read_stats_5p$read_length +
    read_stats_3p$read_length +
    inner_distance

  # -------------------------------
  # Define fragment status
  # -------------------------------
  fstats <- get_fragment_mutation_statuses(mstat_1 = read_info_1$mstat, mstat_2 = read_info_2$mstat)

  # -------------------------------
  # Build an adaptative dataframe
  # -------------------------------
  final_row_fragment <- list(
    Chromosome             = chr,
    Position               = pos,
    Ref                    = ref,
    Alt                    = alt,
    Input_Mutation         = input_mutation_info,
    Fragment_Id            = fragment_name,
    Fragment_QC            = "OK",
    Fragment_Status_Simple = fstats$Simple,
    Fragment_Status_Detail = fstats$Detail,
    Fragment_Size          = absolute_size,
    Inner_Distance         = inner_distance,
    Read_5p                = identity_5p,
    Read_5p_Status         = read_info_5p$mstat,
    Read_3p_Status         = read_info_3p$mstat,
    MAPQ_5p                = read_stats_5p$MAPQ,
    MAPQ_3p                = read_stats_3p$MAPQ,
    BASE_5p                = read_info_5p$base,
    BASE_3p                = read_info_3p$base,
    BASQ_5p                = read_info_5p$basq,
    BASQ_3p                = read_info_3p$basq,
    CIGAR_5p               = read_stats_5p$CIGAR,
    CIGAR_3p               = read_stats_3p$CIGAR,
    POS_5p                 = read_stats_5p$POS,
    POS_3p                 = read_stats_3p$POS
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

    final_row_fragment$Nb_Fragment_Bases_Softclip_5p <-
      fragment_bases_softclip_5p_3p$nb_softclip_5p
    final_row_fragment$Nb_Fragment_Bases_Softclip_3p <-
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
