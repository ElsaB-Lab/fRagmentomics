# Project : ElsaBLab_fRagmentomics

#' @title Identify read1 and read2 from flags
#' @description Determines which sequencing read is the first mate (read1) and which is the second mate (read2).
#'
#' @param df_fragment_reads A dataframe with two rows, each representing a sequencing read.
#' @return A named list containing two elements: `read1` (first mate) and `read2` (second mate).
#' 
#' @noRd
identify_reads <- function(df_fragment_reads) {
  # Extract flag values
  flag_a <- bitvalues_from_bam_flag(
    as.integer(df_fragment_reads[1, "flag"]),
    bitnames = c("isFirstMateRead","isSecondMateRead")
  )
  flag_b <- bitvalues_from_bam_flag(
    as.integer(df_fragment_reads[2, "flag"]),
    bitnames = c("isFirstMateRead","isSecondMateRead")
  )
  
  if ((flag_a[, "isFirstMateRead"] == 1) & (flag_b[, "isSecondMateRead"] == 1)) {
    list(read1 = df_fragment_reads[1, ], read2 = df_fragment_reads[2, ])
  } else if ((flag_b[, "isFirstMateRead"] == 1) &
             (flag_a[, "isSecondMateRead"] == 1)) {
    list(read1 = df_fragment_reads[2, ], read2 = df_fragment_reads[1, ])
  } else {
    stop("Invalid read flags. One read should be first mate, the other second mate.")
  }
}

#' @title Extract read-level statistics
#' @description Retrieves key attributes from a single sequencing read.
#'
#' @param df_read A dataframe containing a single read (one row).
#' @return A list with extracted read properties: mapq, tlen, cigar, pos, query, qual, and read_length.
#' 
#' @noRd
get_read_stats <- function(df_read) {
  list(
    mapq        = df_read$mapq,
    tlen        = df_read$tlen,
    cigar       = df_read$cigar,
    pos         = df_read$pos,
    query       = df_read$seq,
    qual        = df_read$qual,
    read_length = nchar(df_read$seq)
  )
}

#' @title Process a single sequencing fragment
#' @description Extracts and processes relevant data from a sequencing fragment (paired-end reads).
#'
#' @param df_sam A dataframe containing sequencing reads.
#' @param fragment_name Name of the fragment (paired-end reads).
#' @param sample_id Sample identifier.
#' @param chr Chromosome of interest.
#' @param pos Genomic position of interest.
#' @param ref Reference base.
#' @param alt Alternative base (for SNV) or sequence (for insertion).
#' @param mutation_type Type of mutation: "mutation", "deletion", or "insertion".
#'
#' @return A dataframe with the processed fragment information.
#' @importFrom stringr str_extract
#' 
#' @noRd 
process_fragment <- function(df_sam,
                             fragment_name,
                             sample,
                             chr,
                             pos,
                             ref,
                             alt,
                             mutation_type,
                             report_bases_fragment_5p_3p) {
  
  df_fragment_reads <- df_sam[df_sam[,1,drop=TRUE] == fragment_name, , drop=FALSE]
  
  # -------------------------------
  # Sanity check fragments
  # -------------------------------
  fragment_qc <- ""
  if (!all(df_fragment_reads[,3] == chr)) {
    fragment_qc <- paste(fragment_qc, "& Read(s) from another chromosome than", chr)
  }
  if (nrow(df_fragment_reads) != 2) {
    fragment_qc <- paste(fragment_qc, "&", nrow(df_fragment_reads), "read(s)")
  }
  
  if (fragment_qc != "") {
    return(data.frame(
      Chromosome     = chr,
      Position       = pos,
      Fragment       = fragment_name,
      Fragment_Check = fragment_qc
    ))
  }
  
  # Separate read1/read2
  reads       <- identify_reads(df_fragment_reads)
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
  inner_distance <- get_fragment_inner_distance(
    read1_stats$pos, read1_stats$cigar, read1_stats$read_length,
    read2_stats$pos, read2_stats$cigar, read2_stats$read_length
  )
  absolute_size <- read1_stats$read_length +
                   read2_stats$read_length +
                   inner_distance
  
  # -------------------------------
  # Define fragment status
  # -------------------------------
  fragment_status <- process_fragment_status(mutation_type, r_info1$base, r_info2$base)

  # -------------------------------
  # Define 3' and 5' reads
  # -------------------------------
  identity_3p <- ifelse(read1_stats$pos <= read2_stats$pos, 2, 1)
  identity_5p <- ifelse(identity_3p == 2, 1, 2)

  if (identity_5p == 1){
    read_stats_5p <- read1_stats
    read_stats_3p <- read2_stats
  } else {
    read_stats_5p <- read2_stats
    read_stats_3p <- read1_stats
  }

  # -------------------------------
  # Build an adaptative dataframe 
  # -------------------------------
  final_columns <- list(
    Chromosome       = chr,
    Position         = pos,
    Ref              = ref,
    Alt              = alt,
    Sample_Id        = sample,
    Fragment_Id      = fragment_name,
    Fragment_QC      = "OK",
    Fragment_Mutated = fragment_status,
    Absolute_size    = absolute_size,
    Inner_distance   = inner_distance,
    Read_5p          = identity_5p,
    MAPQ_5p          = read_stats_5p$mapq,
    MAPQ_3p          = read_stats_3p$mapq,
    ALT_5p           = read_stats_5p$base,
    ALT_3p           = read_stats_3p$base,
    BASQ_5p          = read_stats_5p$qual,
    BASQ_3p          = read_stats_3p$qual,
    TLEN             = abs(read1_stats$tlen)
  )

  # -------------------------------
  # Define n base 5' and 3'
  # -------------------------------
  # if report_5p_3p_bases_fragment != 0, we call the function process_get_fragment_bases
  if (report_5p_3p_bases_fragment != 0) {
    process_get_fragment_bases_5p_3p <- process_get_fragment_bases(
      read_stats_5p$query, 
      read_stats_3p$query,
      read_stats_5p$qual,
      read_stats_3p$qual
    )
    
    final_columns$Fragment_bases_5p  <- process_get_fragment_bases_5p_3p$Fragment_bases_5p
    final_columns$Fragment_bases_3p  <- process_get_fragment_bases_5p_3p$Fragment_bases_3p
    final_columns$Fragment_Qbases_5p <- process_get_fragment_bases_5p_3p$Fragment_Qbases_5p
    final_columns$Fragment_Qbases_3p <- process_get_fragment_bases_5p_3p$Fragment_Qbases_3p
  }

  # -------------------------------
  # Define number of soft clipped bases in 5' and 3'
  # -------------------------------
  # If report_softclip is TRUE, we call the function process_get_fragment_bases_softclip
  if (report_softclip) {
    process_get_fragment_bases_softclip_5p_3p <- process_get_fragment_bases_softclip(
      read_stats_5p$cigar,
      read_stats_3p$cigar
    )
    
    final_columns$Nb_fragment_bases_softclip_5p <- process_get_fragment_bases_softclip_5p_3p$Nb_SC_5p
    final_columns$Nb_fragment_bases_softclip_3p <- process_get_fragment_bases_softclip_5p_3p$Nb_SC_3p
  }

  # Creation of the dataframe with the wanted columns
  result_df <- as.data.frame(final_columns)

  return(result_df)
}
