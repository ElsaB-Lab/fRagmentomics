skip_if_no_bcftools()

test_that("apply_bcftools_norm", {
  # Create synthetic fasta for chr17
  # Context constructed (TG/TGT/TGTAGG)
  prefix <- paste(rep("A", 125), collapse = "")
  # Sequence at pos 126: T, 127: G, 128: T, 129: A, 130: G, 131: G
  roi <- "TGTAGG"
  suffix <- paste(rep("A", 50), collapse = "")
  seq_chr17 <- paste0(prefix, roi, suffix)

  fasta_env <- setup_test_fasta(c("chr17" = seq_chr17))
  on.exit(cleanup_test_fasta(fasta_env))

  # Valid cases
  expect_equal(
    apply_bcftools_norm("chr17", 126, "TG", "TGT", fasta_env$path, tempdir(), verbose = TRUE),
    data.frame(chr = "chr17", pos = 127, ref = "G", alt = "GT", stringsAsFactors = FALSE),
    ignore_attr = TRUE
  )

  expect_equal(
    apply_bcftools_norm("chr17", 126, "TGTA", "TG", fasta_env$path, tempdir(), verbose = TRUE),
    data.frame(chr = "chr17", pos = 127, ref = "GTA", alt = "G", stringsAsFactors = FALSE),
    ignore_attr = TRUE
  )

  expect_equal(
    apply_bcftools_norm("chr17", 126, "TGT", "TGTAGG", fasta_env$path, tempdir(), verbose = TRUE),
    data.frame(chr = "chr17", pos = 128, ref = "T", alt = "TAGG", stringsAsFactors = FALSE),
    ignore_attr = TRUE
  )

  expect_equal(
    apply_bcftools_norm("chr17", 126, "TG", "TG", fasta_env$path, tempdir(), verbose = TRUE),
    data.frame(chr = "chr17", pos = 126, ref = "TG", alt = "TG", stringsAsFactors = FALSE),
    ignore_attr = TRUE
  )

  # Error cases
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "AGT", "A", fasta_env$path, tempdir(), verbose = TRUE)))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "TAT", "TA", fasta_env$path, tempdir(), verbose = TRUE)))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "T", "", fasta_env$path, tempdir(), verbose = TRUE)))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "T", "-", fasta_env$path, tempdir(), verbose = TRUE)))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "", "GG", fasta_env$path, tempdir(), verbose = TRUE)))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, ".", "G", fasta_env$path, tempdir(), verbose = TRUE)))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "-", "GG", fasta_env$path, tempdir(), verbose = TRUE)))
})

test_that("apply_bcftools_norm is quiet when verbose = FALSE", {
  seq_chr17 <- paste0(paste(rep("A", 125), collapse = ""), "TGTAGG")
  fasta_env <- setup_test_fasta(c("chr17" = seq_chr17))
  on.exit(cleanup_test_fasta(fasta_env))

  msgs <- testthat::capture_messages(
    res <- apply_bcftools_norm("chr17", 126, "TG", "TG", fasta_env$path, tempdir(), verbose = FALSE)
  )

  expect_equal(length(msgs), 0, info = "No messages should be emitted when verbose = FALSE")

  expect_equal(
    res,
    data.frame(chr = "chr17", pos = 126, ref = "TG", alt = "TG", stringsAsFactors = FALSE),
    ignore_attr = TRUE
  )
})
