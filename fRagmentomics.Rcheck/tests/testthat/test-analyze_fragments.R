# A helper to set up file paths consistently across tests
setup_test_files <- function() {
  list(
    mut = system.file("testdata/mutations/", "mutations_cfdna-test-01_chr1_27433000_27435000.tsv", package = "fRagmentomics"),
    bam = system.file("testdata/bam/", "cfdna-test-01_chr1_27433000_27435000.bam", package = "fRagmentomics"),
    fasta = system.file("testdata/fasta/hg19/", "hg19_chr1_27433000_27435000.fa", package = "fRagmentomics")
  )
}

test_that("analyze_fragments works on standard inputs and returns a data frame", {
  files <- setup_test_files()
  temp_out <- tempdir()

  # Execute the function
  results <- analyze_fragments(
    mut = files$mut,
    bam = files$bam,
    fasta = files$fasta,
    sample_id = "test_sample",
    output_folder = temp_out,
    n_cores = 2
  )

  # --- 1. Check the R object output ---
  expect_s3_class(results, "data.frame")
  expect_gt(nrow(results), 0)

  expected_cols <- c("Chromosome", "Position", "Ref", "Alt", "Fragment_Id", "VAF")
  expect_true(all(expected_cols %in% colnames(results)))

  # --- 2. Check file output ---
  expected_file <- file.path(temp_out, "test_sample_df_fRagmentomics.tsv")
  expect_true(file.exists(expected_file))

  # Clean up
  file.remove(expected_file)
})

test_that("analyze_fragments handles different parameter settings correctly", {
  files <- setup_test_files()

  # Test with additional reporting columns enabled
  results_extra_cols <- analyze_fragments(
    mut = files$mut,
    bam = files$bam,
    fasta = files$fasta,
    report_tlen = TRUE,
    report_softclip = TRUE,
    n_cores = 1
  )

  # TLEN is already in the default output, but we confirm others
  extra_cols <- c("Nb_Fragment_Bases_Softclip_5p", "Nb_Fragment_Bases_Softclip_3p")
  expect_true(all(extra_cols %in% colnames(results_extra_cols)))

  # Test with string-based mutation input
  results_string_mut <- analyze_fragments(
    mut = "chr1:685:C:T",
    bam = files$bam,
    fasta = files$fasta,
    n_cores = 1
  )
  expect_s3_class(results_string_mut, "data.frame")
  expect_gt(nrow(results_string_mut), 0)
  expect_equal(results_string_mut$Position[1], 685)
})


test_that("analyze_fragments handles file output logic correctly", {
  files <- setup_test_files()
  temp_out <- tempdir()

  # --- 1. No sample_id provided ---
  default_filename <- "df_fRagmentomics.tsv"
  expected_file_default <- file.path(temp_out, default_filename)
  if (file.exists(expected_file_default)) file.remove(expected_file_default)

  # Expect a message about the default filename being used
  expect_message(
    analyze_fragments(
      mut = files$mut, bam = files$bam, fasta = files$fasta,
      sample_id = NA, output_folder = temp_out, n_cores = 1L
    ),
    "No sample_id provided, using default filename"
  )
  expect_true(file.exists(expected_file_default))

  # --- 2. File overwrite message ---
  expect_message(
    analyze_fragments(
      mut = files$mut, bam = files$bam, fasta = files$fasta,
      sample_id = NA, output_folder = temp_out, n_cores = 1L
    ),
    "already exists and will be overwritten"
  )

  # Clean up
  file.remove(expected_file_default)

  # --- 3. No output_folder provided ---
  # Ensure no file is written
  results_no_output <- analyze_fragments(
    mut = files$mut, bam = files$bam, fasta = files$fasta,
    sample_id = "no_output_test", output_folder = NA, n_cores = 1L
  )
  expect_false(file.exists(file.path(temp_out, "no_output_test_df_fRagmentomics.tsv")))
  expect_s3_class(results_no_output, "data.frame") # Still returns an R object
})


test_that("analyze_fragments handles edge cases gracefully", {
  # Add this skip to avoid errors if bcftools isn't found
  files <- setup_test_files()

  # --- 1. Mutation with no reads should cause the function to stop ---
  mut_no_reads <- "chr1:1:A:T"

  # Test for the specific error message
  suppressWarnings({
    expect_error(
      analyze_fragments(
        mut = mut_no_reads,
        bam = files$bam,
        fasta = files$fasta,
        n_cores = 1
      ),
      "The final dataframe post fRagmentomics is empty."
    )
  })

  # --- 2. Mix of valid and no-read mutations should give a warning and return valid results ---
  mut_mixed_df <- data.frame(
    CHROM = c("chr1", "chr1"),
    POS = c(685, 1),
    REF = c("C", "A"),
    ALT = c("T", "T")
  )
  mut_mixed_file <- tempfile(fileext = ".tsv")
  write.table(mut_mixed_df, mut_mixed_file, sep = "\t", row.names = FALSE, quote = FALSE)
  on.exit(unlink(mut_mixed_file))

  # This call should produce a warning but NOT an error
  expect_warning(
    results_mixed <- analyze_fragments(
      mut = mut_mixed_file,
      bam = files$bam,
      fasta = files$fasta,
      n_cores = 1
    ),
    "No read covers the position of interest"
  )

  # The output should contain results ONLY for the valid mutation
  expect_s3_class(results_mixed, "data.frame")
  expect_gt(nrow(results_mixed), 0)
  expect_true(all(results_mixed$Position == 685))
})
