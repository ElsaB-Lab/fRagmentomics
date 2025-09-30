#' @title Process a single DNA fragment from paired-end reads
#'
#' @description This is a high-level worker function that orchestrates the complete analysis of a single DNA fragment
#' (represented by a pair of reads) at a specific variant locus. It performs quality control, identifies reads, extracts
#' features, determines mutation status, and returns all results in a structured format.
#'
#' @details
#' The function executes the following pipeline for each fragment:
#' 1.  It subsets the reads from 'df_sam' that match the given 'fragment_name'.
#' 2.  It performs initial quality control checks (e.g., proper pairing,
#'     chromosome consistency) via 'process_fragment_reads_qc'.
#' 3.  It identifies the 5' (forward strand) and 3' (reverse strand) reads
#'     based on their SAM FLAGs.
#' 4.  If 'remove_softclip' is 'TRUE', it trims soft-clipped bases from the
#'     sequences, qualities, and CIGAR strings.
#' 5.  It calls 'get_base_basq_mstat_from_read' to determine the allele and
#'     mutation status for each individual read at the variant position.
#' 6.  It calculates the precise fragment size using 'get_fragment_size'.
#' 7.  It consolidates the two read statuses into a final fragment status
#'     (e.g., "MUT", "DISCORDANT") using 'get_mutation_status_of_fragment'.
#' 8.  It assembles and returns a single-row data frame containing all
#'     extracted information.
#'
#' @inheritParams run_fRagmentomics
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
#' @importFrom Rsamtools bamFlagAsBitMatrix
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
                                      remove_softclip,
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
  fragment_qc <- process_fragment_reads_qc(df_fragment_reads, chr)

  # If the fragment fails QC, return a dataframe with NAs
  if (fragment_qc != "") {
    result_df <- create_empty_fragment_row(
      chr, pos, ref, alt, input_mutation_info, fragment_name, fragment_qc,
      sample_id, report_tlen, report_5p_3p_bases_fragment, report_softclip
    )
    return(result_df)
  }

  # -------------------------------
  # Define 3' and 5' reads
  # -------------------------------
  # Get a numeric matrix of FLAG attributes for both reads in the fragment.
  flag_matrix <- bamFlagAsBitMatrix(df_fragment_reads$FLAG)

  # Sanity check that the fragment is a valid pair.
  # It must contain exactly one "first mate" read.
  if (sum(flag_matrix[, "isFirstMateRead"]) != 1) {
    stop(sprintf("Fragment '%s' is not a valid R1/R2 pair.", fragment_name))
  }
  # It must contain one forward and one reverse read.
  if (sum(flag_matrix[, "isMinusStrand"]) != 1) {
    stop(sprintf("Fragment '%s' does not have one forward and one reverse read.", fragment_name))
  }

  # Identify the row index of the 5p read (forward strand, where isMinusStrand is 0).
  idx_5p <- which(flag_matrix[, "isMinusStrand"] == 0)

  # The 3p read is the other one (reverse strand, where isMinusStrand is 1).
  idx_3p <- which(flag_matrix[, "isMinusStrand"] == 1)

  # Get read bam info for read 5' and read 3'
  read_stats_5p <- get_read_stats(df_fragment_reads[idx_5p, ])
  read_stats_3p <- get_read_stats(df_fragment_reads[idx_3p, ])


  # -------------------------------
  # If remove_softclip = TRUE, remove softclip for analysis
  # -------------------------------
  if (remove_softclip) {
    read_5p_info_without_softclip <- remove_softclip(read_stats_5p)
    read_3p_info_without_softclip <- remove_softclip(read_stats_3p)

    # Check if trimming resulted in an empty sequence or CIGAR
    if (read_5p_info_without_softclip$SEQ == "" || read_5p_info_without_softclip$CIGAR == "") {
      # The read is invalid after trimming, so fail the whole fragment.
      result_df <- create_empty_fragment_row(
        chr, pos, ref, alt, input_mutation_info, fragment_name,
        fragment_qc = "Read empty after softclip trim",
        sample_id, report_tlen, report_5p_3p_bases_fragment, report_softclip
      )
      return(result_df)
    }

    read_stats_5p$SEQ <- read_5p_info_without_softclip$SEQ
    read_stats_5p$QUAL <- read_5p_info_without_softclip$QUAL
    read_stats_5p$CIGAR <- read_5p_info_without_softclip$CIGAR
    read_stats_5p$read_length <- read_5p_info_without_softclip$read_length
  }

  # -------------------------------
  # Get read sequence, read base qualities, and read mutation status
  # -------------------------------
  read_info_5p <- get_base_basq_mstat_from_read(
    chr, pos, ref, alt, read_stats_5p, fasta_fafile, fasta_seq
  )
  read_info_3p <- get_base_basq_mstat_from_read(
    chr, pos, ref, alt, read_stats_3p, fasta_fafile, fasta_seq
  )

  # -------------------------------
  # Compute fragment size
  # -------------------------------
  fragment_size <- get_fragment_size(read_stats_5p, read_stats_3p)

  # -------------------------------
  # Define fragment status
  # -------------------------------
  fstats <- get_mutation_status_of_fragment(mstat_5p = read_info_5p$mstat, mstat_3p = read_info_3p$mstat)

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
    Fragment_Status_Simple = as.character(fstats$Simple),
    Fragment_Status_Detail = as.character(fstats$Detail),
    Fragment_Size          = as.integer(fragment_size),
    Read_5p_Status         = as.character(read_info_5p$mstat),
    Read_3p_Status         = as.character(read_info_3p$mstat),
    FLAG_5p                = as.integer(read_stats_5p$FLAG),
    FLAG_3p                = as.integer(read_stats_3p$FLAG),
    MAPQ_5p                = as.integer(read_stats_5p$MAPQ),
    MAPQ_3p                = as.integer(read_stats_3p$MAPQ),
    BASE_5p                = as.character(read_info_5p$base),
    BASE_3p                = as.character(read_info_3p$base),
    BASQ_5p                = as.character(read_info_5p$basq),
    BASQ_3p                = as.character(read_info_3p$basq),
    CIGAR_5p               = as.character(read_stats_5p$CIGAR),
    CIGAR_3p               = as.character(read_stats_3p$CIGAR),
    POS_5p                 = as.integer(read_stats_5p$POS),
    POS_3p                 = as.integer(read_stats_3p$POS)
  )

  # -------------------------------
  # Put sample if not NA
  # -------------------------------
  if (!is.na(sample_id)) {
    as.character(final_row_fragment$Sample_Id <- sample_id)
  }

  # -------------------------------
  # Put sample if not NA
  # -------------------------------
  if (report_tlen) {
    if (is.null(read_stats_5p$TLEN) || is.na(read_stats_5p$TLEN)) {
      warning("TLEN is NULL or missing, cannot be reported.")
      final_row_fragment$TLEN <- "Warning: TLEN is NULL"
    } else {
      as.integer(final_row_fragment$TLEN <- abs(read_stats_5p$TLEN))
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

    as.character(final_row_fragment$Fragment_Bases_5p <-
      fragment_bases_5p_3p$fragment_bases_5p)
    as.character(final_row_fragment$Fragment_Bases_3p <-
      fragment_bases_5p_3p$fragment_bases_3p)
    as.character(final_row_fragment$Fragment_Basqs_5p <-
      fragment_bases_5p_3p$fragment_basqs_5p)
    as.character(final_row_fragment$Fragment_Basqs_3p <-
      fragment_bases_5p_3p$fragment_basqs_3p)
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

    as.integer(final_row_fragment$Nb_Fragment_Bases_Softclip_5p <-
      fragment_bases_softclip_5p_3p$nb_softclip_5p)
    as.integer(final_row_fragment$Nb_Fragment_Bases_Softclip_3p <-
      fragment_bases_softclip_5p_3p$nb_softclip_3p)
  }

  # Add VAF column
  final_row_fragment$VAF <- NA_real_

  # Creation of the dataframe with the wanted columns
  result_df <- as.data.frame(final_row_fragment)
  result_df
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
    FLAG        = df_read$FLAG,
    MAPQ        = df_read$MAPQ,
    TLEN        = df_read$TLEN,
    CIGAR       = df_read$CIGAR,
    POS         = df_read$POS,
    SEQ         = df_read$SEQ,
    QUAL        = df_read$QUAL,
    read_length = nchar(df_read$SEQ)
  )
}
