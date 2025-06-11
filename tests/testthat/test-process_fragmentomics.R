test_that("Process fRagmentomics", {
  mut_test <- system.file("testdata/mutations/", "mutations_cfdna-test-01_chr1_27433000_27434000.tsv", package = "fRagmentomics")
  bam_test <- system.file("testdata/bam/", "cfdna-test-01_chr1_27433000_27434000.bam", package = "fRagmentomics")
  fasta_test <- system.file("testdata/fasta/hg19/", "hg19_chr1_27433000_27434000.fa", package = "fRagmentomics")

  # Read the mutation file
  mutation_file_expected <- read.table(
    mut_test,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  # Execute the function with sample_id = NA
  res1 <- process_fragmentomics(
    mut = mut_test,
    bam = bam_test,
    fasta = fasta_test,
    sample_id = NA,
    neg_offset_mate_search = -200,
    pos_offset_mate_search = 200,
    one_based = TRUE,
    flag_keep = 0x03,
    flag_remove = 0x900,
    report_tlen = TRUE,
    report_softclip = TRUE,
    report_5p_3p_bases_fragment = 5,
    tmp_folder = tempdir(),
    output_file = "./test.tsv", # file will be created in the current directory
    n_cores = 2
  )

  # Check that each row of res1 matches a line in the mutation file
  for (i in seq_len(nrow(res1))) {
    expect_true(
      any(
        mutation_file_expected$CHROM == res1$Chromosome[i] &
          mutation_file_expected$POS == res1$Position[i] &
          mutation_file_expected$REF == res1$Ref[i] &
          mutation_file_expected$ALT == res1$Alt[i]
      )
    )
  }

  # Check that we have the expected columns
  expected_cols <- c(
    "Chromosome", "Position", "Ref", "Alt",
    "Fragment_QC", "Fragment_Status", "Fragment_size",
    "Inner_distance", "MAPQ_5p", "MAPQ_3p", "Read_5p",
    "BASE_5p", "BASE_3p", "BASQ_5p", "BASQ_3p",
    "Pos_bam_3p", "Pos_bam_5p", "CIGAR_3p", "CIGAR_5p",
    "VAF"
  )

  expect_true(
    all(expected_cols %in% colnames(res1)),
    info = "Some expected columns are missing in the output dataframe."
  )

  # Check that the file exists
  output_file <- "./test.tsv"
  expect_true(
    file.exists(output_file),
    info = "Output file was not created."
  )

  # Try to remove it
  file_removed <- file.remove(output_file)
  expect_true(
    file_removed,
    info = paste("Failed to remove output file", output_file)
  )
  expect_false(
    file.exists(output_file),
    info = paste("File still exists after removal attempt:", output_file)
  )
})
