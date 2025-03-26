test_that("apply_bcftools_norm", {
  fasta_38 <- system.file("testdata/fasta/hg38", "hg38_chr17_7676001_7676400.fa", package = "fRagmentomics")

  # Both fastas are in "chrX" format
  # Load fasta as FaFile
  fasta_loaded_38 <- FaFile(fasta_38)
  open(fasta_loaded_38)

  # Valid cases
  expect_equal(
    apply_bcftools_norm("chr17", 126, "TG", "TGT", fasta_38),
    data.frame(chr = "chr17", pos = 127, ref = "G", alt = "GT", stringsAsFactors = FALSE),
    ignore_attr = TRUE
  )

  expect_equal(
    apply_bcftools_norm("chr17", 126, "TGTA", "TG", fasta_38),
    data.frame(chr = "chr17", pos = 127, ref = "GTA", alt = "G", stringsAsFactors = FALSE),
    ignore_attr = TRUE
  )

  expect_equal(
    apply_bcftools_norm("chr17", 126, "TGT", "TGTAGG", fasta_38),
    data.frame(chr = "chr17", pos = 128, ref = "T", alt = "TAGG", stringsAsFactors = FALSE),
    ignore_attr = TRUE
  )

  # expect_equal(
  #     apply_bcftools_norm("chr17", 126, "TGTA", "AG", fasta_38),
  #     data.frame(
  #     chr = c("chr17", "chr17"),
  #     pos = c(126, 127),
  #     ref = c("T", "GTA"),
  #     alt = c("A", "G"),
  #     stringsAsFactors = FALSE
  #     ),
  #     ignore_attr = TRUE
  # )

  # expect_equal(
  #     apply_bcftools_norm("chr17", 125, "GTGT", "GTATCC", fasta_38),
  #     data.frame(
  #     chr = c("chr17", "chr17"),
  #     pos = c(127, 128),
  #     ref = c("G", "T"),
  #     alt = c("A", "TCC"),
  #     stringsAsFactors = FALSE
  #     ),
  #     ignore_attr = TRUE
  # )

  expect_equal(
    apply_bcftools_norm("chr17", 126, "TG", "TG", fasta_38),
    data.frame(chr = "chr17", pos = 126, ref = "TG", alt = "TG", stringsAsFactors = FALSE),
    ignore_attr = TRUE
  )

  # Error cases
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "AGT", "A", fasta_38)))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "TAT", "TA", fasta_38)))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "T", "", fasta_38)))
  # expect_null(apply_bcftools_norm("chr17", 126, "T", ".", fasta_38))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "T", "-", fasta_38)))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "", "GG", fasta_38)))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, ".", "G", fasta_38)))
  expect_null(suppressWarnings(apply_bcftools_norm("chr17", 126, "-", "GG", fasta_38)))
})
