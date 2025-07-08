#' @title Analyze fragments
#'
#' @param mut Path to a .vcf or .tsv file or string representation chr:pos:ref:alt of a mutation.
#' @param bam Path to a BAM file.
#' @param fasta Path to the FASTA file for the reference sequence used for generating the BAM file.
#' @param sample_id Sample identifier.
#' @param neg_offset_mate_search Integer. Use in read_bam.
#'  Represents the number of nucleotides to extend upstream (negative direction)
#'  from the position of interest when querying the BAM file with Rsamtools.
#'  his extension ensures that paired reads are retrieved, even if only one mate
#'  overlaps the queried position.
#' @param pos_offset_mate_search Integer. Use in read_bam.
#' @param one_based Boolean. TRUE if fasta is in one based. False if in 0 based.
#' @param flag_bam_list A named list of logicals for filtering reads based on their SAM flag
#'   NA = Filter is ignored, TRUE = The read MUST have this flag, FALSE = The read MUST NOT have this flag.
#' @param report_tlen Boolean. Whether to include the TLEN (template length)
#'  information in the output.
#' @param report_softclip Boolean. Whether to include the number of soft-clipped
#'  bases at the fragment extremities in the output.
#' @param report_5p_3p_bases_fragment Integer. Whether to include N fragment
#'  extremity bases in the output.
#' @param cigar_free_indel_match Boolean. If activated, the information from the CIGAR is disregarded when determining the
#'  mutation status of a read for indel. Instead the mutation status is determined by comparing the sequence of the read
#'  to the sequence of the wild-type reference and the mutated reference. Activating this option may lead to discordant
#'  genotyping of reads compared to the information provided by the CIGAR for indels. On the other hand, when
#'  activated, it may rescue mutated genotypes for indel that would be missed in cases where the representation of the
#'  indel in the CIGAR does not match the norm of bcftools of the mutation being analyzed.
#' @param tmp_folder Character vector for the folder temporary path.
#' @param output_file Character vector for the output file path.
#' @param n_cores Number of cores for parallel computation.
#'
#' @return A dataframe containing extracted fragment-level information.
#'
#' @importFrom Rsamtools FaFile
#' @importFrom utils write.table
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @importFrom progressr with_progress progressor
#'
#' @export
analyze_fragments <- function(
    mut,
    bam,
    fasta,
    sample_id = NA,
    neg_offset_mate_search = -1000,
    pos_offset_mate_search = 1000,
    one_based = TRUE,
    flag_bam_list = list(
      isPaired = TRUE,
      isProperPair = NA,
      isUnmappedQuery = FALSE,
      hasUnmappedMate = NA,
      isMinusStrand = NA,
      isMateMinusStrand = NA,
      isFirstMateRead = NA,
      isSecondMateRead = NA,
      isSecondaryAlignment = FALSE,
      isSupplementaryAlignment = FALSE,
      isNotPassingQualityControls = NA,
      isDuplicate = FALSE
    ),
    report_tlen = FALSE,
    report_softclip = FALSE,
    report_5p_3p_bases_fragment = 5,
    cigar_free_indel_match = FALSE,
    tmp_folder = tempdir(),
    output_file = NA,
    n_cores = 8) {
  # Load inputs, check parameters and normalize ========================================================================

  # Ensure expected parameters are treated as integers (and not double)
  neg_offset_mate_search <- as.integer(neg_offset_mate_search)
  pos_offset_mate_search <- as.integer(pos_offset_mate_search)
  report_5p_3p_bases_fragment <- as.integer(report_5p_3p_bases_fragment)
  n_cores <- as.integer(n_cores)

  # Check if bam and fasta exist
  # Check if fasta is indexed
  check_parameters(
    mut,
    bam,
    fasta,
    sample_id,
    neg_offset_mate_search,
    pos_offset_mate_search,
    one_based,
    flag_bam_list,
    report_tlen,
    report_softclip,
    report_5p_3p_bases_fragment,
    cigar_free_indel_match,
    tmp_folder,
    output_file,
    n_cores
  )

  # Load and remove bad mutations
  df_mut_raw <- read_mut(mut)
  df_mut_raw <- remove_bad_mut(df_mut_raw)

  # Load fasta as FaFile
  fasta_fafile <- Rsamtools::FaFile(fasta)
  open(fasta_fafile)
  # Close fasta_fafile at the end of the function
  on.exit(close(fasta_fafile), add = TRUE)

  # Normalize mutations
  df_mut_norm <- normalize_mut(df_mut_raw, fasta, fasta_fafile, one_based, tmp_folder)

  # Run per-mutation analysis ==========================================================================================
  # Initialize parallel cluster
  future::plan(future::multisession, workers = n_cores)
  # on.exit() close the parallelisation at the end
  on.exit(future::plan("sequential"), add = TRUE)

  # Create final df
  df_fragments_info_final <- data.frame()

  # Loop on each row of the mut_info
  for (i in seq_len(nrow(df_mut_norm))) {
    chr_norm <- df_mut_norm[i, "chr"]
    pos_norm <- df_mut_norm[i, "pos"]
    ref_norm <- df_mut_norm[i, "ref"]
    alt_norm <- df_mut_norm[i, "alt"]
    input_mutation_info <- df_mut_norm[i, "input_mutation_info"]

    # Read and extract bam around the mutation position
    # Return a truncated sam
    df_sam <- read_bam(
      bam,
      chr_norm,
      pos_norm,
      neg_offset_mate_search,
      pos_offset_mate_search,
      flag_bam_list
    )

    # Check read_bam
    if (is.null(df_sam)) {
      warning_message <- paste0(
        "No read covers the position of interest for the mutation ",
        chr_norm, ":", pos_norm, ":", ref_norm, ">", alt_norm, ". Skipping."
      )
      warning(warning_message, immediate. = TRUE, call. = FALSE)
      next
    }

    # Make one large request to FASTA to avoid doing time-consuming requests to FASTA for each fragment
    # Calculations below serve to ensure we extract just enough reference sequence for all the per-fragment
    # requests
    if (nchar(ref_norm) == nchar(alt_norm)) {
      motif_len <- nchar(alt_norm)
    } else {
      motif_len <- max(nchar(ref_norm) - 1, nchar(alt_norm) - 1)
    }

    max_read_seq_len <- max(nchar(df_sam$SEQ))
    min_read_pos <- min(df_sam[df_sam$RNAME == chr_norm, "POS"])
    max_read_pos <- max(df_sam[df_sam$RNAME == chr_norm, "POS"])
    max_fetch_len_ref <- max_read_seq_len + motif_len

    # -1 for the max value of n_match_base_before used later in the code
    min_fetch_pos <- min_read_pos - 1
    # +1 is the max value of n_match_base_after used later in the code
    max_fetch_pos <- max_read_pos + max_read_seq_len + motif_len + 1

    # get sequence from fasta
    fetch_seq <- get_seq_from_fasta(chr_norm, min_fetch_pos, max_fetch_pos, fasta_fafile)
    fasta_seq <- list(chr = chr_norm, start = min_fetch_pos, end = max_fetch_pos, seq = fetch_seq)

    # Process fragmentomics on truncated bam and the mutation
    # Extract unique fragment names
    fragments_names <- unique(df_sam[, 1, drop = TRUE])
    n_fragments <- length(fragments_names)

    progressr::with_progress({
      # Creation of a progressor to show the progression
      p <- progressr::progressor(steps = n_fragments)

      # future_lapply is the equivalent of lapply in parallele (manage the necessary export)
      results_list <- future.apply::future_lapply(
        fragments_names,
        function(fragment_name) {
          # Tell progressor a step is done
          p()
          extract_fragment_features(
            df_sam                      = df_sam,
            fragment_name               = fragment_name,
            sample_id                   = sample_id,
            chr                         = chr_norm,
            pos                         = pos_norm,
            ref                         = ref_norm,
            alt                         = alt_norm,
            report_tlen                 = report_tlen,
            report_softclip             = report_softclip,
            report_5p_3p_bases_fragment = report_5p_3p_bases_fragment,
            cigar_free_indel_match      = cigar_free_indel_match,
            fasta_seq                   = fasta_seq,
            input_mutation_info         = input_mutation_info
          )
        }
      )
    })

    # future_lapply return a list. Combine in one dataframe.
    df_fragments_info <- do.call(rbind, results_list)

    # Calculate VAF of the fragment
    if (any(df_fragments_info$Fragment_Status_Simple == "MUT", na.rm = TRUE)) {
      total_mut <- sum(df_fragments_info$Fragment_Status_Simple == "MUT", na.rm = TRUE)
      total_non_target_mut <- sum(df_fragments_info$Fragment_Status_Simple == "NON-TARGET MUT", na.rm = TRUE)

      if (total_mut + total_non_target_mut == 0) {
        df_fragments_info$VAF <- 0
      } else {
        df_fragments_info$VAF <- 100 * total_mut / (total_mut + total_non_target_mut)
      }
    } else {
      if (all(is.na(df_fragments_info$Fragment_Status_Simple == "MUT"), na.rm = TRUE)) {
        df_fragments_info$VAF <- NA
      } else {
        df_fragments_info$VAF <- 0
      }
    }

    # Fusion into the final df
    df_fragments_info_final <- rbind(df_fragments_info_final, df_fragments_info)
  }

  # Check if the df post fRagmentomics is not empty
  if (nrow(df_fragments_info_final) == 0) {
    stop("The final dataframe post fRagmentomics is empty.")
  }

  # -------------------------------
  # Write output file if output directory specified
  # -------------------------------
  if (!(is.na(output_file) || output_file == "")) {
    # Write file
    write.table(
      df_fragments_info_final,
      output_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  # Return final dataframe
  df_fragments_info_final
}
