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
                             mutation_type) {
  
  df_fragment_reads <- df_sam[df_sam[,1,drop=TRUE] == fragment_name, , drop=FALSE]
  
  # Sanity checks
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
  
  # Retrieve mutation information
  r_info1 <- get_mutation_info(mutation_type, pos, ref, alt, read1_stats)
  r_info2 <- get_mutation_info(mutation_type, pos, ref, alt, read2_stats)
  
  # Compute distances and fragment size
  n_unique_del_ins <- get_nb_unique_del_ins(
    read1_stats$pos, read1_stats$cigar,
    read2_stats$pos, read2_stats$cigar
  )
  inner_distances <- get_fragment_inner_distance(
    read1_stats$pos, read1_stats$cigar, read1_stats$read_length,
    read2_stats$pos, read2_stats$cigar, read2_stats$read_length
  )
  absolute_size <- read1_stats$read_length +
                   read2_stats$read_length +
                   inner_distances$b

  identity_RIGHT <- ifelse(read1_stats$pos <= read2_stats$pos, 2, 1)
  identity_LEFT  <- ifelse(identity_RIGHT == 2, 1, 2)
  
  return(data.frame(
    Sample_Id        = sample_id,
    Chromosome       = chr,
    Position         = pos,
    Ref              = ref,
    Alt              = alt,
    Fragment         = fragment_name,
    Fragment_Check   = "OK",
    Mutation_type    = mutation_type,
    Mapq_Read_1      = read1_stats$mapq,
    Mapq_Read_2      = read2_stats$mapq,
    Base_Read_1      = r_info1$base,
    Base_Read_2      = r_info2$base,
    Qual_Read_1      = r_info1$qual,
    Qual_Read_2      = r_info2$qual,
    TLEN             = abs(read1_stats$tlen),
    Absolute_size    = absolute_size
  ))
}
