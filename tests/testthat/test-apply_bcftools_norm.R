# This helper will skip all tests in a file if bcftools is missing

skip_if_no_bcftools <- function() {
  if (Sys.which("bcftools") == "") {
    testthat::skip("bcftools not available for testing")
  }
}

skip_if_no_bcftools()

test_that("apply_bcftools_norm", {
  fasta_38 <- system.file("testdata/fasta/hg38", "hg38_chr17_7676001_7676400.fa", package = "fRagmentomics")

  # Both fastas are in "chrX" format
  # Load fasta as FaFile
  fasta_loaded_38 <- FaFile(fasta_38)
  open(fasta_loaded_38)

  # Valid cases
  expect_equal(
    apply_bcftools_norm("chr17", 126, "TG", "TGT", fasta_38, tempdir(), verbose = TRUE),
    data.frame(chr = "chr17", pos = 127, ref = "G", alt = "GT", stringsAsFactors = FALSE),
    ignore_attr = TRUE
  )

  expect_equal(
    apply_bcftools_norm("chr17", 126, "TGTA", "TG", fasta_38, tempdir(), verbose = TRUE),
    data.frame(chr = "chr17", pos = 127, ref = "GTA", alt = "G", stringsAsFactors = FALSE),
    ignore_attr = TRUE
  )

  expect_equal(
    apply_bcftools_norm("chr17", 126, "TGT", "TGTAGG", fasta_38, tempdir(), verbose = TRUE),
    data.frame(chr = "chr17", pos = 128, ref = "T", alt = "TAGG", stringsAsFactors = FALSE),
    ignore_attr = TRUE
  )

  expect_equal(
    apply_bcftools_norm("chr17", 126, "TG", "TG", fasta_38, tempdir(), verbose = TRUE),
    data.frame(chr = "chr17", pos = 126, ref = "TG", alt = "TG", stringsAsFactors = FALSE),
    ignore_attr = TRUE
  )

  # Error cases
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "AGT", "A", fasta_38, tempdir(), verbose = TRUE)))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "TAT", "TA", fasta_38, tempdir(), verbose = TRUE)))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "T", "", fasta_38, tempdir(), verbose = TRUE)))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "T", "-", fasta_38, tempdir(), verbose = TRUE)))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "", "GG", fasta_38, tempdir(), verbose = TRUE)))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, ".", "G", fasta_38, tempdir(), verbose = TRUE)))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "-", "GG", fasta_38, tempdir(), verbose = TRUE)))

  close(fasta_loaded_38)
})

test_that("apply_bcftools_norm is quiet when verbose = FALSE", {
  msgs <- testthat::capture_messages(
    res <- apply_bcftools_norm("chr17", 126, "TG", "TG", fasta_38, tempdir(), verbose = FALSE
    )
  )

  expect_equal(length(msgs), 0, info = "No messages should be emitted when verbose = FALSE")

  expect_equal(
    res,
    data.frame(chr = "chr17", pos = 126, ref = "TG", alt = "TG", stringsAsFactors = FALSE),
    ignore_attr = TRUE
  )
})
