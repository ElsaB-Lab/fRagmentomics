# tests/test-run_fRagmentomics.R

# Helper to locate package testdata reliably across tests
setup_test_files <- function() {
  list(
    mut   = system.file("testdata/mutations/", "mutations_cfdna-test-01_chr17_7576000_7579000.tsv", package = "fRagmentomics"),
    bam   = system.file("testdata/bam/", "cfdna-test-01_chr17_7576000_7579000.bam", package = "fRagmentomics"),
    fasta = system.file("testdata/fasta/hg19/", "hg19_chr17_7576000_7579000.fa", package = "fRagmentomics")
  )
}

test_that("returns NULL when writing to file and creates parent directories if needed", {
  files <- setup_test_files()

  # Build a nested, non-existing directory path under tempdir()
  nested_dir <- file.path(tempdir(), "fragmentomics_test", "deep", "subdir")
  out_file <- file.path(nested_dir, sprintf("fragmentomics_%s.tsv", as.integer(runif(1, 1, 1e6))))

  # Ensure the parent directory does not exist beforehand
  if (dir.exists(nested_dir)) unlink(nested_dir, recursive = TRUE, force = TRUE)
  expect_false(dir.exists(nested_dir))

  # Capture all messages in one run and assert on their content
  msgs <- testthat::capture_messages(
    res <- run_fRagmentomics(
      mut = files$mut, bam = files$bam, fasta = files$fasta,
      sample_id = "test_sample",
      output_path = out_file,
      verbose = TRUE,
      n_cores = 4L
    )
  )

  # Check both message kinds and return value (use any() over grepl on the vector)
  expect_true(any(grepl("Folder '.+' does not exist and will be created", msgs)))
  expect_true(any(grepl("^Writing results to:", msgs)))
  expect_null(res)
  expect_true(file.exists(out_file))

  # Second call should emit an overwrite message (capture again)
  msgs2 <- testthat::capture_messages(
    res2 <- run_fRagmentomics(
      mut = files$mut,
      bam = files$bam,
      fasta = files$fasta,
      output_path = out_file,
      verbose = TRUE,
      n_cores = 4L
    )
  )
  expect_true(any(grepl("already exists and will be overwritten", msgs2)))
  expect_true(any(grepl("^Writing results to:", msgs2)))
  expect_null(res2)

  # Clean up
  unlink(file.path(tempdir(), "fragmentomics_test"), recursive = TRUE, force = TRUE)
})

test_that("is quiet in verbose = FALSE (no stray messages) and still writes output", {
  files <- setup_test_files()

  # Fresh nested path under tempdir()
  nested_dir <- file.path(tempdir(), "fragmentomics_test_quiet", "deep", "subdir")
  out_file <- file.path(nested_dir, sprintf("fragmentomics_%s.tsv", as.integer(runif(1, 1, 1e6))))

  if (dir.exists(nested_dir)) unlink(nested_dir, recursive = TRUE, force = TRUE)
  expect_false(dir.exists(nested_dir))

  # Capture only messages; warnings (if any) are not captured here.
  msgs <- testthat::capture_messages(
    res <- run_fRagmentomics(
      mut = files$mut, bam = files$bam, fasta = files$fasta,
      sample_id = "test_sample",
      output_path = out_file,
      verbose = FALSE, # quiet mode
      n_cores = 4L
    )
  )

  # Assert no messages were emitted in quiet mode
  expect_equal(length(msgs), 0, info = "No messages should be emitted when verbose = FALSE")

  # Function still returns NULL and writes the file
  expect_null(res)
  expect_true(file.exists(out_file))

  # Second call (overwrite) should also remain silent in quiet mode
  msgs2 <- testthat::capture_messages(
    res2 <- run_fRagmentomics(
      mut = files$mut, bam = files$bam, fasta = files$fasta,
      sample_id = "test_sample",
      output_path = out_file,
      verbose = FALSE,
      n_cores = 4L
    )
  )
  expect_equal(length(msgs2), 0, info = "No messages on overwrite when verbose = FALSE")
  expect_null(res2)

  # Clean up
  unlink(file.path(tempdir(), "fragmentomics_test_quiet"), recursive = TRUE, force = TRUE)
})

test_that("returns a data frame when output_path is NA / NULL / '' and writes nothing", {
  files <- setup_test_files()

  # 1) output_path = NA
  res_na <- run_fRagmentomics(
    mut = files$mut,
    bam = files$bam,
    fasta = files$fasta,
    sample_id = "no_output_na",
    output_path = NA,
    verbose = TRUE,
    n_cores = 4L
  )
  expect_s3_class(res_na, "data.frame")
  expect_false(any(grepl("no_output_na", list.files(tempdir(), full.names = TRUE))))

  # 2) output_path = NULL
  res_null <- run_fRagmentomics(
    mut = files$mut,
    bam = files$bam,
    fasta = files$fasta,
    sample_id = "no_output_null",
    output_path = NULL,
    verbose = TRUE,
    n_cores = 4L
  )
  expect_s3_class(res_null, "data.frame")
  expect_false(any(grepl("no_output_null", list.files(tempdir(), full.names = TRUE))))

  # 3) output_path = "" (empty string)
  res_empty <- run_fRagmentomics(
    mut = files$mut,
    bam = files$bam,
    fasta = files$fasta,
    sample_id = "no_output_empty",
    output_path = "",
    verbose = TRUE,
    n_cores = 4L
  )
  expect_s3_class(res_empty, "data.frame")
  expect_false(any(grepl("no_output_empty", list.files(tempdir(), full.names = TRUE))))
})

test_that("adds expected columns when reporting flags are enabled and supports string mutation input", {
  files <- setup_test_files()

  # With extra reporting columns enabled
  results_extra <- run_fRagmentomics(
    mut = files$mut,
    bam = files$bam,
    fasta = files$fasta,
    report_tlen = TRUE,
    report_softclip = TRUE,
    verbose = TRUE,
    n_cores = 4L
  )
  expect_s3_class(results_extra, "data.frame")
  expect_true(all(c("Nb_Fragment_Bases_Softclip_5p", "Nb_Fragment_Bases_Softclip_3p") %in% colnames(results_extra)))

  # String-based mutation input
  results_str <- run_fRagmentomics(
    mut = "chr17:2191:T:C",
    bam = files$bam,
    fasta = files$fasta,
    verbose = TRUE,
    n_cores = 4L
  )
  expect_s3_class(results_str, "data.frame")
  expect_gt(nrow(results_str), 0)
  expect_true("Position" %in% names(results_str))
  expect_true(any(results_str$Position == 2191))
})

test_that("emits QC removal message unless retain_fail_qc = TRUE", {
  files <- setup_test_files()

  # Default (retain_fail_qc = FALSE): expect removal message among messages
  msgs_rm <- testthat::capture_messages(
    run_fRagmentomics(
      mut = files$mut,
      bam = files$bam,
      fasta = files$fasta,
      verbose = TRUE,
      n_cores = 4L
    )
  )
  expect_true(any(grepl("Removed .* fragments that fail quality checks", msgs_rm)))

  # With retain_fail_qc = TRUE: the removal message must be absent
  msgs_keep <- testthat::capture_messages(
    run_fRagmentomics(
      mut = files$mut,
      bam = files$bam,
      fasta = files$fasta,
      retain_fail_qc = TRUE,
      verbose = TRUE,
      n_cores = 4L
    )
  )
  expect_false(any(grepl("Removed .* fragments that fail quality checks", msgs_keep)))
})

test_that("gracefully handles loci with no covering reads and mixed-valid inputs", {
  files <- setup_test_files()

  # 1) Single mutation with no reads -> final DF empty -> error
  mut_no_reads <- "chr17:1:A:T"
  suppressWarnings({
    expect_error(
      run_fRagmentomics(
        mut = mut_no_reads,
        bam = files$bam,
        fasta = files$fasta,
        verbose = TRUE,
        n_cores = 4L
      ),
      regexp = "^The final fRagmentomic dataframe is empty\\.$"
    )
  })

  # 2) Mixed mutations: one valid (685) + one with no reads (1).
  mut_mixed_df <- data.frame(
    CHROM = c("chr17", "chr17"),
    POS   = c(2191, 1),
    REF   = c("T", "T"),
    ALT   = c("C", "A")
  )
  mut_mixed_file <- tempfile(fileext = ".tsv")
  write.table(mut_mixed_df, mut_mixed_file, sep = "\t", row.names = FALSE, quote = FALSE)
  on.exit(unlink(mut_mixed_file), add = TRUE)

  # Capture warning text while keeping the actual return value
  captured_warnings <- character(0)

  res_mixed <- withCallingHandlers(
    run_fRagmentomics(
      mut = mut_mixed_file,
      bam = files$bam,
      fasta = files$fasta,
      verbose = TRUE,
      n_cores = 1L
    ),
    warning = function(w) {
      captured_warnings <<- c(captured_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("No read covers the position of interest", captured_warnings, fixed = TRUE)))

  expect_s3_class(res_mixed, "data.frame")
  expect_gt(nrow(res_mixed), 0)
  expect_true(all(res_mixed$Position == 2191))
})

test_that("VAF is computed from status labels when present", {
  files <- setup_test_files()

  res <- run_fRagmentomics(
    mut = files$mut,
    bam = files$bam,
    fasta = files$fasta,
    verbose = TRUE,
    n_cores = 4L
  )
  expect_true("VAF" %in% names(res))

  # VAF should be numeric in [0, 100] or NA when all statuses are NA
  vaf <- res$VAF
  expect_true(is.numeric(vaf))
  expect_true(all(is.na(vaf) | (vaf >= 0 & vaf <= 100)))
})
