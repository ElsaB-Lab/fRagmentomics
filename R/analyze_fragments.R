#' Setup parallel computations
#'
#' @param n_cores Number of cores to use for parallel processing.
#'
#' @return A parallel cluster object.
#'
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#'
#' @noRd
setup_parallel_computations <- function(n_cores) {
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  cl
}


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
#' @param flag_keep Character vector. Use in read_bam.
#'  Represents the SAM flags that should be kept when filtering alignments.
#' @param flag_remove Character vector. Use in read_bam.
#'  Represents the SAM flags that should be excluded when filtering alignments.
#' @param report_tlen Boolean. Whether to include the TLEN (template length)
#'  information in the output.
#' @param report_softclip Boolean. Whether to include the number of soft-clipped
#'  bases at the fragment extremities in the output.
#' @param report_5p_3p_bases_fragment Integer. Whether to include N fragment
#'  extremity bases in the output.
#' @param cigar_free_mode Boolean. If activated, the information from the CIGAR is disregarded when determining the
#'  mutation status of a read. Instead the mutation status is determined by comparing the sequence of the read to the
#'  sequence of the wild-type reference and the mutated reference. Activating this option may lead to discordant
#'  genotyping of reads compared to the information provided by the CIGAR for indels. On the other hand, when
#'  activated, it may rescue mutated genotypes for indel that would be missed in cases where the representation of the
#'  indel in the CIGAR does not match the norm of bcftools of the mutation being analyzed.
#' @param tmp_folder Character vector for the folder temporary path.
#' @param output_file Character vector for the output file path. Mandatory.
#' @param n_cores Number of cores for parallel computation.
#'
#' @return A dataframe containing extracted fragment-level information.
#'
#' @importFrom Rsamtools FaFile
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#' @importFrom parallel stopCluster
#' @importFrom utils write.table
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
    flag_keep = 0x03,
    flag_remove = 0x900,
    report_tlen = FALSE,
    report_softclip = FALSE,
    report_5p_3p_bases_fragment = 5,
    cigar_free_mode = FALSE,
    tmp_folder = tempdir(),
    output_file = NA,
    n_cores = 1) {
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
    one_based, flag_keep,
    flag_remove,
    report_tlen,
    report_softclip,
    report_5p_3p_bases_fragment,
    cigar_free_mode,
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

  # Normalize mutations
  df_mut_norm <- normalize_mut(df_mut_raw, fasta, fasta_fafile, one_based, tmp_folder)

  # Run per-mutation analysis ==========================================================================================
  # Initialize parallel cluster
  cl <- setup_parallel_computations(n_cores)

  # Create final df
  df_fragments_info_final <- data.frame()

  # Loop on each row of the mut_info
  for (i in seq_len(nrow(df_mut_norm))) {
    chr_norm <- df_mut_norm[i, "chr"]
    pos_norm <- df_mut_norm[i, "pos"]
    ref_norm <- df_mut_norm[i, "ref"]
    alt_norm <- df_mut_norm[i, "alt"]

    # Read and extract bam around the mutation position
    # Return a truncated sam
    df_sam <- read_bam(
      bam,
      chr_norm,
      pos_norm,
      neg_offset_mate_search,
      pos_offset_mate_search,
      flag_keep,
      flag_remove
    )

    # Process fragmentomics on truncated bam and the mutation
    # Extract unique fragment names
    fragments_names <- unique(df_sam[, 1, drop = TRUE])
    n_fragments <- length(fragments_names)

    j <- NULL # to avoid "no visible binding for global variable" in R CMD check
    # Parallel execution
    df_fragments_info <- foreach::foreach(
      j = 1:n_fragments,
      .combine = rbind,
      .inorder = FALSE,
      .multicombine = FALSE,
      .packages = "fRagmentomics"
    ) %do% {
      extract_fragment_features(
        df_sam                      = df_sam,
        fragment_name               = fragments_names[j],
        sample_id                   = sample_id,
        chr                         = chr_norm,
        pos                         = pos_norm,
        ref                         = ref_norm,
        alt                         = alt_norm,
        report_tlen                 = report_tlen,
        report_softclip             = report_softclip,
        report_5p_3p_bases_fragment = report_5p_3p_bases_fragment,
        cigar_free_mode             = cigar_free_mode,
        fasta_fafile                = fasta_fafile
      )
    }

    # -------------------------------
    # Calculate VAF of the fragment
    # -------------------------------
    df_fragments_info$VAF <- NA
    if (any(df_fragments_info$Fragment_Status == "MUT")) {
      total_mut <- sum(df_fragments_info$Fragment_Status == "MUT", na.rm = TRUE)
      total <- sum(df_fragments_info$Fragment_Status == "Non-target MUT", na.rm = TRUE)

      df_fragments_info$VAF <- 100 * total_mut / total
    } else {
      df_fragments_info$VAF <- NA
    }

    # Fusion into the final df
    df_fragments_info_final <- rbind(df_fragments_info_final, df_fragments_info)
  }

  # Stop cluster
  parallel::stopCluster(cl)

  # Close fasta
  close(fasta_fafile)

  # Check if the df post fRagmentomics is not empty
  if (nrow(df_fragments_info_final) == 0) {
    stop("The final dataframe post fRagmentomics is empty.")
  }

  # -------------------------------
  # Write output file
  # -------------------------------
  # Write file
  write.table(
    df_fragments_info_final,
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # Return final dataframe
  df_fragments_info_final
}
