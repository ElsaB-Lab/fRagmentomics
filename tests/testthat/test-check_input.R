test_that("check_input", {
  # ---------------------------------------------------------------------------
  # Create temporary test directory and file paths
  # ---------------------------------------------------------------------------
  test_dir <- file.path(tempdir(), "test_data_check_input")
  if (!dir.exists(test_dir)) dir.create(test_dir, recursive = TRUE)

  # Define file paths
  mut_missing <- file.path(test_dir, "missing.vcf")
  mut_existing <- file.path(test_dir, "existing.vcf.gz")
  bam_existing <- file.path(test_dir, "existing.bam")
  bam_missing <- file.path(test_dir, "missing.bam")
  bam_no_index <- file.path(test_dir, "bam_no_index.bam")
  fasta_missing <- file.path(test_dir, "missing.fa")
  fasta_no_index <- file.path(test_dir, "no_index.fa")
  fasta_with_index <- file.path(test_dir, "with_index.fa")

  # ---------------------------------------------------------------------------
  # Prepare fake FASTA files
  # ---------------------------------------------------------------------------
  # Create two FASTA files: one without index and one with an index.
  writeLines(">chr1\nAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT", fasta_no_index)
  writeLines(">chr1\nAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT", fasta_with_index)

  # Index the FASTA that should have an index.
  Rsamtools::indexFa(fasta_with_index)

  # ---------------------------------------------------------------------------
  # Prepare a SAM file and convert it to BAM
  # ---------------------------------------------------------------------------
  sam_file <- file.path(tempdir(), "test.sam")
  writeLines(c(
    "@HD\tVN:1.0\tSO:unsorted",
    "@SQ\tSN:chr1\tLN:1000",
    "read1\t0\tchr1\t1\t255\t10M\t*\t0\t0\tAGCTAGCTAG\t*"
  ), sam_file)

  # Convert SAM to BAM using asBam (creates sorted BAM)
  sorted_bam <- file.path(test_dir, "sorted")
  Rsamtools::asBam(sam_file, sorted_bam, overwrite = TRUE)
  bam_existing <- paste0(sorted_bam, ".bam")

  # Create a BAM without an index by copying the existing BAM.
  bam_no_index <- file.path(test_dir, "bam_no_index.bam")
  file.copy(bam_existing, bam_no_index, overwrite = TRUE)

  # Ensure bam_existing has an index.
  if (!file.exists(paste0(bam_existing, ".bai"))) {
    Rsamtools::indexBam(bam_existing)
  }

  # ---------------------------------------------------------------------------
  # Create an existing mutation file.
  # ---------------------------------------------------------------------------
  file.create(mut_existing)

  # ---------------------------------------------------------------------------
  # Define valid default parameters for the non-file arguments.
  # ---------------------------------------------------------------------------
  valid_params <- list(
    sample                       = "test_sample",
    neg_offset_mate_search       = -1000L,
    pos_offset_mate_search       = 1000L,
    one_based                    = TRUE,
    flag_keep                    = 0x03,
    flag_remove                  = 0x900,
    report_tlen                  = FALSE,
    report_softclip              = FALSE,
    report_5p_3p_bases_fragment  = 5,
    tmp_folder                   = tempdir(),
    n_cores                      = 1
  )

  # ---------------------------------------------------------------------------
  # 1. Valid input: all files exist and required indexes are present.
  # ---------------------------------------------------------------------------
  expect_silent(
    check_input(
      mut = mut_existing,
      bam = bam_existing,
      fasta = fasta_with_index,
      sample = valid_params$sample,
      neg_offset_mate_search = valid_params$neg_offset_mate_search,
      pos_offset_mate_search = valid_params$pos_offset_mate_search,
      one_based = valid_params$one_based,
      flag_keep = valid_params$flag_keep,
      flag_remove = valid_params$flag_remove,
      report_tlen = valid_params$report_tlen,
      report_softclip = valid_params$report_softclip,
      report_5p_3p_bases_fragment = valid_params$report_5p_3p_bases_fragment,
      tmp_folder = valid_params$tmp_folder,
      n_cores = valid_params$n_cores
    )
  )

  # ---------------------------------------------------------------------------
  # 2. Missing mutation file: expect an error.
  # ---------------------------------------------------------------------------
  expect_error(
    check_input(
      mut = mut_missing,
      bam = bam_existing,
      fasta = fasta_with_index,
      sample = valid_params$sample,
      neg_offset_mate_search = valid_params$neg_offset_mate_search,
      pos_offset_mate_search = valid_params$pos_offset_mate_search,
      one_based = valid_params$one_based,
      flag_keep = valid_params$flag_keep,
      flag_remove = valid_params$flag_remove,
      report_tlen = valid_params$report_tlen,
      report_softclip = valid_params$report_softclip,
      report_5p_3p_bases_fragment = valid_params$report_5p_3p_bases_fragment,
      tmp_folder = valid_params$tmp_folder,
      n_cores = valid_params$n_cores
    ),
    "Error: The Mutation file does not exist",
    fixed = TRUE
  )

  # ---------------------------------------------------------------------------
  # 3. Missing BAM file: expect an error.
  # ---------------------------------------------------------------------------
  expect_error(
    check_input(
      mut = mut_existing,
      bam = bam_missing,
      fasta = fasta_with_index,
      sample = valid_params$sample,
      neg_offset_mate_search = valid_params$neg_offset_mate_search,
      pos_offset_mate_search = valid_params$pos_offset_mate_search,
      one_based = valid_params$one_based,
      flag_keep = valid_params$flag_keep,
      flag_remove = valid_params$flag_remove,
      report_tlen = valid_params$report_tlen,
      report_softclip = valid_params$report_softclip,
      report_5p_3p_bases_fragment = valid_params$report_5p_3p_bases_fragment,
      tmp_folder = valid_params$tmp_folder,
      n_cores = valid_params$n_cores
    ),
    "Error: The BAM file does not exist",
    fixed = TRUE
  )

  # ---------------------------------------------------------------------------
  # 4. Missing FASTA file: expect an error.
  # ---------------------------------------------------------------------------
  expect_error(
    check_input(
      mut = mut_existing,
      bam = bam_existing,
      fasta = fasta_missing,
      sample = valid_params$sample,
      neg_offset_mate_search = valid_params$neg_offset_mate_search,
      pos_offset_mate_search = valid_params$pos_offset_mate_search,
      one_based = valid_params$one_based,
      flag_keep = valid_params$flag_keep,
      flag_remove = valid_params$flag_remove,
      report_tlen = valid_params$report_tlen,
      report_softclip = valid_params$report_softclip,
      report_5p_3p_bases_fragment = valid_params$report_5p_3p_bases_fragment,
      tmp_folder = valid_params$tmp_folder,
      n_cores = valid_params$n_cores
    ),
    "Error: The FASTA file does not exist",
    fixed = TRUE
  )

  # ---------------------------------------------------------------------------
  # 5. Missing BAM index: remove the index from bam_no_index and expect a message.
  # ---------------------------------------------------------------------------
  bai_path <- paste0(bam_no_index, ".bai")
  if (file.exists(bai_path)) file.remove(bai_path)

  expect_message(
    check_input(
      mut = mut_existing,
      bam = bam_no_index,
      fasta = fasta_with_index,
      sample = valid_params$sample,
      neg_offset_mate_search = valid_params$neg_offset_mate_search,
      pos_offset_mate_search = valid_params$pos_offset_mate_search,
      one_based = valid_params$one_based,
      flag_keep = valid_params$flag_keep,
      flag_remove = valid_params$flag_remove,
      report_tlen = valid_params$report_tlen,
      report_softclip = valid_params$report_softclip,
      report_5p_3p_bases_fragment = valid_params$report_5p_3p_bases_fragment,
      tmp_folder = valid_params$tmp_folder,
      n_cores = valid_params$n_cores
    ),
    "Creating BAM index..."
  )

  # ---------------------------------------------------------------------------
  # 6. Missing FASTA index: remove the FASTA index from fasta_no_index and expect a message.
  # ---------------------------------------------------------------------------
  fai_path <- paste0(fasta_no_index, ".fai")
  if (file.exists(fai_path)) file.remove(fai_path)

  expect_message(
    check_input(
      mut = mut_existing,
      bam = bam_existing,
      fasta = fasta_no_index,
      sample = valid_params$sample,
      neg_offset_mate_search = valid_params$neg_offset_mate_search,
      pos_offset_mate_search = valid_params$pos_offset_mate_search,
      one_based = valid_params$one_based,
      flag_keep = valid_params$flag_keep,
      flag_remove = valid_params$flag_remove,
      report_tlen = valid_params$report_tlen,
      report_softclip = valid_params$report_softclip,
      report_5p_3p_bases_fragment = valid_params$report_5p_3p_bases_fragment,
      tmp_folder = valid_params$tmp_folder,
      n_cores = valid_params$n_cores
    ),
    "Creating FASTA index..."
  )

  # ---------------------------------------------------------------------------
  # Clean up temporary test files and folders
  # ---------------------------------------------------------------------------
  unlink(test_dir, recursive = TRUE)
  unlink(sam_file)
})
