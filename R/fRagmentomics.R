# Project : ElsaBLab_fRagmentomics

#' Setup parallel computations
#'
#' @inheritParams fRagmentomics
#' @param n_cores Number of cores to use for parallel processing.
#' @return A parallel cluster object.
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#'
#' @noRd
setup_parallel_computations <- function(n_cores) {
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  return(cl)
}

#' @title fRagmentomics
#'
#' @param bam A dataframe containing sequencing reads.
#' @param mut Could be a .vcf, .tsv or chr:pos:ref:alt.
#' @param fasta A reference genome in FASTA format.
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
#' @param n_cores Number of cores for parallel computation.
#'
#' @return A dataframe containing extracted fragment-level information.
#'
#' @importFrom Rsamtools FaFile
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel stopCluster
#'
#' @export
fRagmentomics <- function(
    bam,
    mut,
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
    n_cores = 1) {
  # -------------------------------
  # Check the inputs and load files
  # -------------------------------
  # Check if bam and fasta exist
  # Check if fasta is indexed
  check_input(
    mut, bam, fasta, sample_id, neg_offset_mate_search, pos_offset_mate_search,
    one_based, flag_keep, flag_remove, report_tlen, report_softclip,
    report_5p_3p_bases_fragment, n_cores
  )

  # Read a vcf ou .tsv file
  # Return a df with all the mutation we want to study
  # Look at the format of mut to know if it is a VCF, mut_file, chr:pos:ref:alt
  if (length(mut) == 1) {
    mut_info <- read_mut(mut)
  } else {
    stop("Error: The parameter 'mut' should be a single value,
    not multiple elements.")
  }

  # Sanity check to validate input transformation
  mut_info_checked <- sanity_check_read_mut(mut_info)

  # Load fasta as FaFile
  fasta_loaded <- Rsamtools::FaFile(fasta)
  open(fasta_loaded)

  # -------------------------------
  # Normalisation of Ref, Alt and Pos
  # -------------------------------
  # Initialisation of the normalized variants
  final_variants <- data.frame()

  for (i in seq_len(nrow(mut_info_checked))) {
    chr <- mut_info_checked[i, 1]
    pos <- mut_info_checked[i, 2]
    ref <- mut_info_checked[i, 3]
    alt <- mut_info_checked[i, 4]

    # Normalization user-provided representation into vcf representation
    mut_info_vcf_normalized <- normalize_user_rep_to_vcf_rep(
      chr = chr,
      pos = pos,
      ref = ref,
      alt = alt,
      fasta_loaded = fasta_loaded,
      one_based = one_based
    )

    # Sanity check to see if ref != fasta
    if (is.null(mut_info_vcf_normalized)) {
      next
    }

    with(mut_info_vcf_normalized, {
      chr_norm <- chr
      pos_norm <- pos
      ref_norm <- ref
      alt_norm <- alt
    })

    # Normalization vcf representation with bcftools norm
    mut_info_bcftools_normalized <- apply_bcftools_norm(
      chr = chr_norm,
      pos = pos_norm,
      ref = ref_norm,
      alt = alt_norm,
      fasta = fasta
    )

    # Sanity check to see if bcftools worked properly
    if (is.null(mut_info_bcftools_normalized)) {
      next
    }

    # Append to the final dataframe
    final_variants <- rbind(final_variants, mut_info_bcftools_normalized)
  }

  # -------------------------------
  # Perform fragment analysis
  # -------------------------------
  # Loop on each row of the mut_info
  for (i in seq_len(nrow(final_variants))) {
    chr_final <- final_variants[i, 1]
    pos_final <- final_variants[i, 2]
    ref_final <- final_variants[i, 3]
    alt_final <- final_variants[i, 4]

    # Return the mutation status in SNV, ins, del, MNP
    mutation_type <- define_mutation_status(ref_final, alt_final)

    # Read and extract bam around the mutation position
    # Return a truncated sam
    df_sam <- read_bam(
      bam,
      chr_final,
      pos_final,
      neg_offset_mate_search,
      pos_offset_mate_search,
      flag_keep,
      flag_remove
    )

    # Process fragmentomics on truncated bam and the mutation
    # Extract unique fragment names
    fragments_names <- unique(df_sam[, 1, drop = TRUE])
    n_fragments <- length(fragments_names)

    # Initialize parallel cluster
    cl <- setup_parallel_computations(n_cores)

    # Parallel execution
    df_fragments_info <- foreach::foreach(
      i = seq_len(n_fragments),
      .combine = rbind,
      .inorder = FALSE,
      .packages = c("stringr")
    ) %dopar% {
      # Call the fragment processing function
      process_fragment(
        df_sam = df_sam,
        fragment_name = fragments_names[i],
        sample_id = sample_id,
        chr = chr_final,
        pos = pos_final,
        ref = ref_final,
        alt = alt_final,
        mutation_type = mutation_type,
        report_tlen = report_tlen,
        report_softclip = report_softclip,
        report_5p_3p_bases_fragment = report_5p_3p_bases_fragment
      )
    }

    # Stop cluster
    parallel::stopCluster(cl)
  }

  return(df_fragments_info)
}
