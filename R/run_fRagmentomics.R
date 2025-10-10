#' @title Analyze fragments
#'
#' @description
#' This is the main function of the package. It provides an end-to-end pipeline for analyzing the allelic state of
#' individual DNA fragments covering specific genomic variants. It takes a list of mutations and an aligned sequencing
#' file (BAM) as input, processes each fragment in parallel, and returns a detailed data frame of results.
#'
#' @details
#' The function executes a multi-step workflow for each variant provided in the
#' 'mut' input:
#' \enumerate{
#'   \item **Input Validation**: All parameters are rigorously checked for correctness (e.g., file existence, data types).
#'     Required file indices ('.bai', '.fai') are created automatically if missing.
#'   \item **Variant Normalization**: The input variants are parsed and normalized into a canonical, left-aligned
#'     representation using a combination of VCF-style indel padding and the external 'bcftools norm' command.
#'   \item **BAM Read Extraction**: For each normalized variant, the function efficiently queries the BAM file to
#'     retrieve all read pairs that cover the genomic locus.
#'   \item **Parallel Fragment Processing**: The core analysis is performed in parallel using the 'future' framework.
#'     Each unique DNA fragment is processed by the 'extract_fragment_features' worker function to determine
#'     its size, quality metrics, and mutation status (e.g., "MUT", "WT", "AMB", "N/I").
#'   \item **VAF Calculation**: After all fragments for a variant are processed, the Variant Allele Frequency (VAF) is calculated.
#'   \item **Output Generation**: Results from all variants are aggregated into a single data frame. If a value for
#'     `output_path` is provided, this data frame is also written to a tab-separated file.
#' }
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
#' @param remove_softclip Boolean. If set to TRUE, trim soft-clipped bases from the 5' end of Read 5p and from the 3' end of Read 3p.
#' @param retain_fail_qc Boolean. If set to TRUE, retain fragments that failed the various quality checks in the output.
#' @param apply_bcftools_norm Boolean. If set to TRUE, apply bcftools norm on each input variant to normalize it.
#'  Require that bcftools command is installed and available in the PATH.
#' @param tmp_folder Character vector for the temporary folder path.
#' @param output_path Character vector for the fragmentomics table output path.
#' @param verbose Boolean. If set to TRUE, print all the warnings and the prints.
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
#' @examples
#' # --- 1. Locate Example Files ---
#' # The package includes small example files to demonstrate its functionality.
#' # We locate them using system.file().
#' mut_file <- system.file(
#'   "extdata", "mutations_cfdna-test-01_chr1_27433000_27435000.tsv",
#'   package = "fRagmentomics"
#' )
#' bam_file <- system.file(
#'   "extdata", "cfdna-test-01_chr1_27433000_27435000.bam",
#'   package = "fRagmentomics"
#' )
#' fasta_file <- system.file(
#'   "extdata", "hg19_chr1_27433000_27435000.fa",
#'   package = "fRagmentomics"
#' )
#'
#' # --- 2. Run the Analysis ---
#' # This single call runs the full analysis pipeline on the example data.
#' # The output file is written to a temporary location to avoid cluttering
#' # the working directory. We use n_cores = 1L for examples.
#' results <- run_fRagmentomics(
#'   mut = mut_file,
#'   bam = bam_file,
#'   fasta = fasta_file,
#'   sample_id = "cfdna-test-01",
#'   n_cores = 1L
#' )
#'
#' # --- 3. View the Results ---
#' # Print the first few rows of the output data frame to see the results.
#' print(head(results))
#'
run_fRagmentomics <- function(
    mut,
    bam,
    fasta,
    sample_id = NA_character_,
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
    remove_softclip = FALSE,
    retain_fail_qc = FALSE,
    apply_bcftools_norm = FALSE,
    tmp_folder = tempdir(),
    output_path = NA_character_,
    verbose = FALSE,
    n_cores = 8) {
  # Load inputs, check parameters and normalize ========================================================================

  # Ensure expected parameters are treated as integers (and not double)
  neg_offset_mate_search <- as.integer(neg_offset_mate_search)
  pos_offset_mate_search <- as.integer(pos_offset_mate_search)
  report_5p_3p_bases_fragment <- as.integer(report_5p_3p_bases_fragment)
  n_cores <- as.integer(n_cores)

  # Check if all params pass QC
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
    remove_softclip,
    retain_fail_qc,
    tmp_folder,
    output_path,
    verbose,
    n_cores
  )

  # Check system dependencey
  if (apply_bcftools_norm){
    bcftools_path <- check_bcftools_is_installed()
    if (verbose){
      print(paste("Found bcftools at:", bcftools_path))
    }
  }

  # Load and remove bad mutations
  df_mut_raw <- read_mut(mut)
  df_mut_raw <- remove_bad_mut(df_mut_raw)

  # Load fasta as FaFile
  fasta_fafile <- Rsamtools::FaFile(fasta)
  open(fasta_fafile)
  # Close fasta_fafile at the end of the function
  on.exit(close(fasta_fafile), add = TRUE)

  # Normalize mutations
  df_mut_norm <- normalize_mut(df_mut_raw, fasta, fasta_fafile, one_based, apply_bcftools_norm, tmp_folder, verbose)

  # Run per-mutation analysis ==========================================================================================
  # Initialize parallel cluster
  if (n_cores == 1L) {
    future::plan("sequential")
  } else {
    future::plan(future::multisession, workers = n_cores)
  }
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
    read_bam_out <- read_bam(
      bam,
      chr_norm,
      pos_norm,
      neg_offset_mate_search,
      pos_offset_mate_search,
      flag_bam_list
    )

    df_sam_badly_oriented <- read_bam_out$badly_oriented
    df_sam <- read_bam_out$well_oriented

    # Check if there are badly-oriented reads and create dataframe with NAs if any are found
    if (!is.null(df_sam_badly_oriented)) {
      df_fragments_info_badly_oriented <- data.frame()

      for (fragment_name in unique(df_sam_badly_oriented$QNAME)) {
        fragment_qc <- "Fragment with badly oriented reads"

        df_fragments_info_badly_oriented <- rbind(
          df_fragments_info_badly_oriented,
          create_empty_fragment_row(
            chr = chr_norm,
            pos = pos_norm,
            ref = ref_norm,
            alt = alt_norm,
            input_mutation_info = input_mutation_info,
            fragment_name = fragment_name,
            fragment_qc = fragment_qc,
            sample_id = sample_id,
            report_tlen = report_tlen,
            report_5p_3p_bases_fragment = report_5p_3p_bases_fragment,
            report_softclip = report_softclip
          )
        )
      }

      # Add to the final df
      df_fragments_info_final <- rbind(df_fragments_info_final, df_fragments_info_badly_oriented)
    }

    # Check that there are well-oriented reads to be processed
    if (is.null(df_sam)) {
      warning(
        sprintf(
          "No read covers the position of interest for the mutation %s:%d:%s>%s. Skipping.",
          chr_norm, pos_norm, ref_norm, alt_norm
        ),
        immediate. = TRUE, call. = FALSE
      )
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
            remove_softclip             = remove_softclip,
            fasta_seq                   = fasta_seq,
            input_mutation_info         = input_mutation_info
          )
        }
      )
    })

    # future_lapply return a list. Combine in one dataframe.
    df_fragments_info <- do.call(rbind, results_list)

    # Calculate VAF of the fragment
    status_simple <- df_fragments_info$Fragment_Status_Simple

    if (all(is.na(status_simple))) {
      df_fragments_info$VAF <- NA_real_
    } else {
      denom <- sum(status_simple %in% c("MUT", "WT", "OTH"), na.rm = TRUE)
      num <- sum(status_simple == "MUT", na.rm = TRUE)
      df_fragments_info$VAF <- if (denom == 0) 0 else 100 * num / denom
    }

    # Fusion into the final df
    df_fragments_info_final <- rbind(df_fragments_info_final, df_fragments_info)
  }

  # Check if the final df is not empty
  if (nrow(df_fragments_info_final) == 0) {
    stop("The final fRagmentomic dataframe is empty.")
  }

  # Remove reads failing QC unless the user explicity wants to retain them
  if (!retain_fail_qc) {
    mask_fail_qc <- df_fragments_info_final[["Fragment_QC"]] != "OK"
    nfrag_fail_qc <- sum(mask_fail_qc)
    df_fragments_info_final <- df_fragments_info_final[!mask_fail_qc, , drop = FALSE]
    if (verbose) {
      message(sprintf("Removed '%d' fragments that fail quality checks from output. Set `retain_fail_qc` to TRUE to retain them.", nfrag_fail_qc))
    }
  }

  # -------------------------------
  # Write output file if output directory specified
  # -------------------------------
  # Only proceed with file writing if an output_path is specified
  if (!is.null(output_path) && !is.na(output_path) && output_path != "") {
    # Check if the file already exists and warn the user if it will be overwritten
    if (file.exists(output_path)) {
      if (verbose) {
        message(sprintf("File '%s' already exists and will be overwritten.", output_path))
      }
    } else {
      # Check if the file parent folder already exists
      output_parent <- dirname(output_path)
      if (!dir.exists(output_parent)) {
        if (verbose) {
          message(sprintf("Folder '%s' does not exist and will be created", output_parent))
        }
        dir.create(output_parent, showWarnings = FALSE, recursive = TRUE)
      }
    }

    # Write the data frame to the file
    if (verbose) {
      message(sprintf("Writing results to: %s", output_path))
    }
    write.table(
      df_fragments_info_final,
      file = output_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    return(NULL)
  }

  return(df_fragments_info_final)
}
