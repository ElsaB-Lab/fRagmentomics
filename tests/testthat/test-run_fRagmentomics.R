# Helper to locate package testdata reliably across tests
setup_test_files <- function() {
  list(
    mut = system.file("extdata/mutation",
      "cfdna-egfr-del_chr7_55241864_55243064_10k.mutations.tsv",
      package = "fRagmentomics"
    ),
    bam = system.file("extdata/bam",
      "cfdna-egfr-del_chr7_55241864_55243064_10k.bam",
      package = "fRagmentomics"
    ),
    fasta = system.file("extdata/fasta",
      "hg19_chr7_55231864_55253064.fa",
      package = "fRagmentomics"
    )
  )
}

test_that("Master Run: Handles File I/O, Messages, QC, and Extra Columns (Verbose=TRUE)", {
  files <- setup_test_files()

  # Setup unique output path
  nested_dir <- file.path(tempdir(), "frag_test_master", "subdir")
  out_file <- file.path(nested_dir, "output.tsv")
  if (dir.exists(nested_dir)) unlink(nested_dir, recursive = TRUE, force = TRUE)

  # --- EXECUTION 1---
  # We activate ALL optional logic
  msgs <- testthat::capture_messages(
    res <- run_fRagmentomics(
      mut = files$mut,
      bam = files$bam,
      fasta = files$fasta,
      sample_id = "test_sample",
      output_path = out_file,
      verbose = TRUE,
      n_cores = 2L,
      report_bam_info = TRUE,
      report_softclip = TRUE
    )
  )

  # --- VERIFICATIONS ---

  # 1. File & Folder Creation logic
  expect_true(any(grepl("Folder '.+' does not exist and will be created", msgs)))
  expect_true(any(grepl("^Writing results to:", msgs)))
  expect_null(res) # Should return NULL when writing to file
  expect_true(file.exists(out_file))

  # 2. QC Logic (Default is retain_fail_qc = FALSE)
  expect_true(any(grepl("Removed .* fragments that fail quality checks", msgs)))

  # 3. Data Integrity (Read back the file)
  res_df <- utils::read.table(out_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Check VAF logic
  expect_true("VAF" %in% names(res_df))
  expect_true(is.numeric(res_df$VAF))

  # Check Extra Columns (report_bam_info/softclip)
  expect_true(all(c("Nb_Fragment_Bases_Softclip_5p", "Nb_Fragment_Bases_Softclip_3p") %in% colnames(res_df)))

  # Clean up
  unlink(file.path(tempdir(), "frag_test_master"), recursive = TRUE, force = TRUE)
})

test_that("Quiet Run & In-Memory Return (Verbose=FALSE, Output=NULL)", {
  files <- setup_test_files()

  # --- EXECUTION 2---
  # 1. Quiet mode (no messages)
  # 2. Output=NULL returns data.frame (covers NA and "" logic too as they share the same 'if' block)
  # 3. String input for mutation (valid coordinate from the test file to ensure it works)

  mut_df <- utils::read.table(files$mut, header = TRUE, stringsAsFactors = FALSE)
  valid_mut_str <- sprintf("%s:%d:%s:%s", mut_df$CHROM[1], mut_df$POS[1], mut_df$REF[1], mut_df$ALT[1])

  msgs <- testthat::capture_messages(
    res <- run_fRagmentomics(
      mut = valid_mut_str,
      bam = files$bam,
      fasta = files$fasta,
      output_path = NULL, # Should return object
      verbose = FALSE, # Should be silent
      n_cores = 2L
    )
  )

  # --- VERIFICATIONS ---
  expect_length(msgs, 0) # Quiet mode check
  expect_s3_class(res, "data.frame") # Return check
  expect_gt(nrow(res), 0)
})

test_that("Error Handling: No Reads Found", {
  files <- setup_test_files()

  # --- EXECUTION 3 ---
  # Use a coordinate clearly outside the BAM range (e.g. 1) to trigger empty result error
  mut_no_reads <- "chr7:1:A:T"

  suppressWarnings({
    expect_error(
      run_fRagmentomics(
        mut = mut_no_reads,
        bam = files$bam,
        fasta = files$fasta,
        verbose = FALSE,
        n_cores = 2L
      ),
      regexp = "^The final fRagmentomic dataframe is empty\\.$"
    )
  })
})

skip_if_no_bcftools()

test_that("Bcftools: User parameter apply_bcftools_norm is correctly used", {
  files <- setup_test_files()

  res <- run_fRagmentomics(
    mut = files$mut,
    bam = files$bam,
    fasta = files$fasta,
    verbose = FALSE,
    apply_bcftools_norm = TRUE,
    n_cores = 2L
  )
  expect_s3_class(res, "data.frame")
})
