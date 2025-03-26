test_that("normalize_user_rep_to_vcf_rep", {
  fasta_19 <- system.file("testdata/fasta/hg19", "hg19_chr1_230709546_230710546.fa", package = "fRagmentomics")
  fasta_38 <- system.file("testdata/fasta/hg38", "hg38_chr1_230709546_230710546.fa", package = "fRagmentomics")

  # Both fastas are in "chrX" format
  # Load fasta as FaFile
  fasta_loaded_19 <- FaFile(fasta_19)
  open(fasta_loaded_19)

  fasta_loaded_38 <- FaFile(fasta_38)
  open(fasta_loaded_38)

  # --------------------------------------------------------------------------
  # Case 1: SNP (one-based)
  # --------------------------------------------------------------------------
  out1 <- normalize_user_rep_to_vcf_rep(
    chr = "chr1",
    pos = 503,
    ref = "G",
    alt = "A",
    fasta = fasta_loaded_19,
    one_based = TRUE
  )
  expect_equal(out1$chr, "chr1")
  expect_equal(out1$pos, 503)
  expect_equal(out1$ref, "G")
  expect_equal(out1$alt, "A")

  # --------------------------------------------------------------------------
  # Case 2: Deletion
  # --------------------------------------------------------------------------
  out2 <- normalize_user_rep_to_vcf_rep(
    chr = "chr1",
    pos = 503,
    ref = "GTT",
    alt = "-",
    fasta = fasta_loaded_19,
    one_based = TRUE
  )
  # For a deletion, we expect an anchor base prepended to both REF and ALT,
  # and pos shifts left by 1.
  expect_equal(out2$chr, "chr1")
  expect_equal(out2$pos, 502)
  expect_equal(out2$ref, "GGTT")
  expect_equal(out2$alt, "G")

  # --------------------------------------------------------------------------
  # Case 3: Insertion
  # --------------------------------------------------------------------------
  out3 <- normalize_user_rep_to_vcf_rep(
    chr = "chr1",
    pos = 501,
    ref = "",
    alt = "AT",
    fasta = fasta_loaded_19,
    one_based = TRUE
  )
  expect_equal(out3$chr, "chr1")
  expect_equal(out3$pos, 501)
  expect_equal(out3$ref, "A")
  expect_equal(out3$alt, "AAT")

  # --------------------------------------------------------------------------
  # Case 4: Deletion with '-' in ALT
  # Input: ref = "ATT", alt = "A--", pos = 3
  # --------------------------------------------------------------------------
  out4 <- normalize_user_rep_to_vcf_rep(
    chr = "chr1",
    pos = 502,
    ref = "GGTT",
    alt = "G---",
    fasta = fasta_loaded_19,
    one_based = TRUE
  )
  expect_equal(out4$chr, "chr1")
  expect_equal(out4$pos, 502)
  expect_equal(out4$ref, "GGTT")
  expect_equal(out4$alt, "G")

  # --------------------------------------------------------------------------
  # Case 5: Complex deletion
  # --------------------------------------------------------------------------
  out5 <- normalize_user_rep_to_vcf_rep(
    chr = "chr1",
    pos = 502,
    ref = "GGTT",
    alt = "A",
    fasta = fasta_loaded_19,
    one_based = TRUE
  )
  expect_equal(out5$chr, "chr1")
  expect_equal(out5$pos, 502)
  expect_equal(out5$ref, "GGTT")
  expect_equal(out5$alt, "A")

  # --------------------------------------------------------------------------
  # Case 6: Complex insertion
  # --------------------------------------------------------------------------
  out6 <- normalize_user_rep_to_vcf_rep(
    chr = "chr1",
    pos = 501,
    ref = "A",
    alt = "GAT",
    fasta = fasta_loaded_19,
    one_based = TRUE
  )
  expect_equal(out6$chr, "chr1")
  expect_equal(out6$pos, 501)
  expect_equal(out6$ref, "A")
  expect_equal(out6$alt, "GAT")

  # --------------------------------------------------------------------------
  # Case 7: SNP, but position is given in 0-based coords
  # --------------------------------------------------------------------------
  out7 <- normalize_user_rep_to_vcf_rep(
    chr = "chr1",
    pos = 502, # 0-based
    ref = "GT",
    alt = "AC",
    fasta = fasta_loaded_19,
    one_based = FALSE
  )
  expect_equal(out7$chr, "chr1")
  expect_equal(out7$pos, 503)
  expect_equal(out7$ref, "GT")
  expect_equal(out7$alt, "AC")

  # --------------------------------------------------------------------------
  # Case 8: SNP normal case but in hg38 position
  # --------------------------------------------------------------------------
  out8 <- normalize_user_rep_to_vcf_rep(
    chr = "chr1",
    pos = 503,
    ref = "A",
    alt = "T",
    fasta = fasta_loaded_38,
    one_based = TRUE
  )
  expect_equal(out8$chr, "chr1")
  expect_equal(out8$pos, 503)
  expect_equal(out8$ref, "A")
  expect_equal(out8$alt, "T")

  # --------------------------------------------------------------------------
  # Case 9: Check chromosome transformation
  # --------------------------------------------------------------------------
  out8 <- normalize_user_rep_to_vcf_rep(
    chr = 1,
    pos = 503,
    ref = "G",
    alt = "A",
    fasta = fasta_loaded_19,
    one_based = TRUE
  )
  expect_equal(out8$chr, "chr1")
  expect_equal(out8$pos, 503)
  expect_equal(out8$ref, "G")
  expect_equal(out8$alt, "A")

  # --------------------------------------------------------------------------
  # Case 10: Check if seq ref != fasta
  # --------------------------------------------------------------------------
  expect_warning(
    out9 <- normalize_user_rep_to_vcf_rep(
      chr = "chr1", # Ensure this matches the FASTA convention
      pos = 503,
      ref = "T", # Incorrect reference allele
      alt = "A",
      fasta = fasta_loaded_19,
      one_based = TRUE
    ),
    "Mismatch found"
  )
  expect_null(out9)
})
