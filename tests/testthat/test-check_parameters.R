test_that("check_parameters", {
  # ---------------------------------------------------------------------------
  # Create temporary test directory and file paths
  # ---------------------------------------------------------------------------
  test_dir <- file.path(tempdir(), "test_data_check_parameters")
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
    sample = "test_sample",
    neg_offset_mate_search = -1000L,
    pos_offset_mate_search = 1000L,
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
    report_5p_3p_bases_fragment = 5L,
    cigar_free_indel_match = FALSE,
    remove_softclip = FALSE,
    tmp_folder = tempdir(),
    output_folder = "./test",
    n_cores = 1L
  )

  # ---------------------------------------------------------------------------
  # 1. Valid input: all files exist and required indexes are present.
  # ---------------------------------------------------------------------------
  expect_silent(
    suppressMessages(
      check_parameters(
        mut = mut_existing,
        bam = bam_existing,
        fasta = fasta_with_index,
        sample = valid_params$sample,
        neg_offset_mate_search = valid_params$neg_offset_mate_search,
        pos_offset_mate_search = valid_params$pos_offset_mate_search,
        one_based = valid_params$one_based,
        flag_bam_list = valid_params$flag_bam_list,
        report_tlen = valid_params$report_tlen,
        report_softclip = valid_params$report_softclip,
        report_5p_3p_bases_fragment = valid_params$report_5p_3p_bases_fragment,
        cigar_free_indel_match = valid_params$cigar_free_indel_match,
        remove_softclip = valid_params$remove_softclip,
        tmp_folder = valid_params$tmp_folder,
        output_folder = valid_params$output_folder,
        n_cores = valid_params$n_cores
      )
    )
  )

  # ---------------------------------------------------------------------------
  # 2. Missing mutation file: expect an error.
  # ---------------------------------------------------------------------------
  expect_error(
    check_parameters(
      mut = mut_missing,
      bam = bam_existing,
      fasta = fasta_with_index,
      sample = valid_params$sample,
      neg_offset_mate_search = valid_params$neg_offset_mate_search,
      pos_offset_mate_search = valid_params$pos_offset_mate_search,
      one_based = valid_params$one_based,
      flag_bam_list = valid_params$flag_bam_list,
      report_tlen = valid_params$report_tlen,
      report_softclip = valid_params$report_softclip,
      report_5p_3p_bases_fragment = valid_params$report_5p_3p_bases_fragment,
      cigar_free_indel_match = valid_params$cigar_free_indel_match,
      remove_softclip = valid_params$remove_softclip,
      tmp_folder = valid_params$tmp_folder,
      output_folder = valid_params$output_folder,
      n_cores = valid_params$n_cores
    ),
    "The Mutation file does not exist",
    fixed = TRUE
  )

  # ---------------------------------------------------------------------------
  # 3. Missing BAM file: expect an error.
  # ---------------------------------------------------------------------------
  expect_error(
    check_parameters(
      mut = mut_existing,
      bam = bam_missing,
      fasta = fasta_with_index,
      sample = valid_params$sample,
      neg_offset_mate_search = valid_params$neg_offset_mate_search,
      pos_offset_mate_search = valid_params$pos_offset_mate_search,
      one_based = valid_params$one_based,
      flag_bam_list = valid_params$flag_bam_list,
      report_tlen = valid_params$report_tlen,
      report_softclip = valid_params$report_softclip,
      report_5p_3p_bases_fragment = valid_params$report_5p_3p_bases_fragment,
      cigar_free_indel_match = valid_params$cigar_free_indel_match,
      remove_softclip = valid_params$remove_softclip,
      tmp_folder = valid_params$tmp_folder,
      output_folder = valid_params$output_folder,
      n_cores = valid_params$n_cores
    ),
    "The BAM file does not exist",
    fixed = TRUE
  )

  # ---------------------------------------------------------------------------
  # 4. Missing FASTA file: expect an error.
  # ---------------------------------------------------------------------------
  expect_error(
    check_parameters(
      mut = mut_existing,
      bam = bam_existing,
      fasta = fasta_missing,
      sample = valid_params$sample,
      neg_offset_mate_search = valid_params$neg_offset_mate_search,
      pos_offset_mate_search = valid_params$pos_offset_mate_search,
      one_based = valid_params$one_based,
      flag_bam_list = valid_params$flag_bam_list,
      report_tlen = valid_params$report_tlen,
      report_softclip = valid_params$report_softclip,
      report_5p_3p_bases_fragment = valid_params$report_5p_3p_bases_fragment,
      cigar_free_indel_match = valid_params$cigar_free_indel_match,
      remove_softclip = valid_params$remove_softclip,
      tmp_folder = valid_params$tmp_folder,
      output_folder = valid_params$output_folder,
      n_cores = valid_params$n_cores
    ),
    "The FASTA file does not exist",
    fixed = TRUE
  )

  # ---------------------------------------------------------------------------
  # 5. Missing BAM index: remove the index from bam_no_index and expect a message.
  # ---------------------------------------------------------------------------
  bai_path <- paste0(bam_no_index, ".bai")
  if (file.exists(bai_path)) file.remove(bai_path)

  expect_message(
    check_parameters(
      mut = mut_existing,
      bam = bam_no_index,
      fasta = fasta_with_index,
      sample = valid_params$sample,
      neg_offset_mate_search = valid_params$neg_offset_mate_search,
      pos_offset_mate_search = valid_params$pos_offset_mate_search,
      one_based = valid_params$one_based,
      flag_bam_list = valid_params$flag_bam_list,
      report_tlen = valid_params$report_tlen,
      report_softclip = valid_params$report_softclip,
      report_5p_3p_bases_fragment = valid_params$report_5p_3p_bases_fragment,
      cigar_free_indel_match = valid_params$cigar_free_indel_match,
      remove_softclip = valid_params$remove_softclip,
      tmp_folder = valid_params$tmp_folder,
      output_folder = valid_params$output_folder,
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
    check_parameters(
      mut = mut_existing,
      bam = bam_existing,
      fasta = fasta_no_index,
      sample = valid_params$sample,
      neg_offset_mate_search = valid_params$neg_offset_mate_search,
      pos_offset_mate_search = valid_params$pos_offset_mate_search,
      one_based = valid_params$one_based,
      flag_bam_list = valid_params$flag_bam_list,
      report_tlen = valid_params$report_tlen,
      report_softclip = valid_params$report_softclip,
      report_5p_3p_bases_fragment = valid_params$report_5p_3p_bases_fragment,
      cigar_free_indel_match = valid_params$cigar_free_indel_match,
      remove_softclip = valid_params$remove_softclip,
      tmp_folder = valid_params$tmp_folder,
      output_folder = valid_params$output_folder,
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

test_that("check_parameters individual parameter validations", {
  # Sample
  expect_error(check_sample(""), "sample ID cannot be empty")
  expect_silent(check_sample(NA))

  # Offsets
  expect_error(check_neg_offset_mate_search(-1000), "must be integer")
  expect_error(check_pos_offset_mate_search(1000), "must be interger")

  expect_error(check_neg_offset_mate_search(as.integer(1000)), "must be negative")
  expect_error(check_pos_offset_mate_search(as.integer(-1000)), "must be positive")

  # Flags
  # A standard, valid partial list
  expect_no_error(
    check_flag_bam_list(list(isPaired = TRUE, isDuplicate = FALSE))
  )
  # A valid list including NA
  expect_no_error(
    check_flag_bam_list(list(isProperPair = NA, isUnmappedQuery = FALSE))
  )
  # An empty list should also be valid
  expect_no_error(
    check_flag_bam_list(list())
  )
  # Test case for non-list input
  expect_error(
    check_flag_bam_list("not a list"),
    "must be a list"
  )
  # Test case for non-logical values in the list
  expect_error(
    check_flag_bam_list(list(isPaired = TRUE, isDuplicate = "FALSE")),
    "must be logical"
  )
  expect_error(
    check_flag_bam_list(list(isUnmappedQuery = 0)),
    "must be logical"
  )
  # Test case for an unnamed list element
  expect_error(
    check_flag_bam_list(list(TRUE, isPaired = FALSE)),
    "must be named"
  )
  # Test case for an invalid flag name
  expect_error(
    check_flag_bam_list(list(isPaired = TRUE, isADuplicate = FALSE)),
    "Invalid name(s) found",
    fixed = TRUE
  )
  # Test case for multiple invalid flag names
  expect_error(
    check_flag_bam_list(list(isPaired = TRUE, isADuplicate = FALSE, isGood = TRUE)),
    "Invalid name(s) found",
    fixed = TRUE
  )

  # Logical flags
  expect_error(check_one_based("TRUE"), "must be a single logical")
  expect_error(check_report_tlen(1L), "must be a single logical")
  expect_error(check_report_softclip(NULL), "must be a single logical")
  expect_error(check_cigar_free_indel_match(NULL), "must be a single logical")
  expect_error(check_remove_softclip(NULL), "must be a single logical")

  # Report 5p/3p bases
  expect_error(check_report_bases_fragm_5p_3p(5), "must be integer")
  expect_error(check_report_bases_fragm_5p_3p(as.integer(-1)), "must be non-negative")

  # Temporary folder
  expect_error(check_tmp_folder(TRUE), "must be a single character")

  # Output folder
  expect_silent(check_output_folder(NA))
  expect_silent(check_output_folder(""))
  expect_error(check_output_folder(42), "'output_folder' must be a single character string.")

  # n_cores
  expect_error(check_n_cores(2), "must be integer")
  expect_error(check_n_cores(as.integer(0)), "must be positive")
})
