test_that("normalize_to_vcf_rep", {
  # ----------------------------------------------------------------------------------------------
  # SETUP: Synthetic Reference
  # Constructed to match the specific logic of the original HG19 snippet at pos 500+
  # Pattern required for tests:
  # Pos 501: A
  # Pos 502: G
  # Pos 503: G
  # Pos 504: T
  # Pos 505: T
  # Sequence segment: ...AGGTT...
  # ----------------------------------------------------------------------------------------------
  seq_chr1 <- paste0(
    strrep("N", 500), # Padding 1-500
    "AGGTT", # ROI 501-505
    strrep("N", 100) # Padding suffix
  )
  names(seq_chr1) <- "chr1"

  fasta_env <- setup_test_fasta(seq_chr1)
  fasta_loaded <- fasta_env$fa_obj

  # --------------------------------------------------------------------------
  # Case 1: SNP (one-based)
  # --------------------------------------------------------------------------
  out1 <- normalize_to_vcf_rep(
    chr = "chr1",
    pos = 503,
    ref = "G",
    alt = "A",
    fasta_fafile = fasta_loaded,
    one_based = TRUE,
    verbose = TRUE
  )
  expect_equal(out1$chr, "chr1")
  expect_equal(out1$pos, 503)
  expect_equal(out1$ref, "G")
  expect_equal(out1$alt, "A")

  # --------------------------------------------------------------------------
  # Case 2: Deletion
  # --------------------------------------------------------------------------
  out2 <- normalize_to_vcf_rep(
    chr = "chr1",
    pos = 503,
    ref = "GTT",
    alt = "-",
    fasta_fafile = fasta_loaded,
    one_based = TRUE,
    verbose = TRUE
  )
  # For a deletion, we expect an anchor base prepended to both REF and ALT,
  # and pos shifts left by 1.
  # Anchor at 502 is 'G'.
  expect_equal(out2$chr, "chr1")
  expect_equal(out2$pos, 502)
  expect_equal(out2$ref, "GGTT")
  expect_equal(out2$alt, "G")

  # --------------------------------------------------------------------------
  # Case 3: Insertion
  # --------------------------------------------------------------------------
  out3 <- normalize_to_vcf_rep(
    chr = "chr1",
    pos = 501,
    ref = "",
    alt = "AT",
    fasta_fafile = fasta_loaded,
    one_based = TRUE,
    verbose = TRUE
  )
  # Anchor at 501 is 'A'.
  expect_equal(out3$chr, "chr1")
  expect_equal(out3$pos, 501)
  expect_equal(out3$ref, "A")
  expect_equal(out3$alt, "AAT")

  # --------------------------------------------------------------------------
  # Case 4: Deletion with '-' in ALT
  # Input: ref = "ATT", alt = "A--", pos = 3
  # --------------------------------------------------------------------------
  out4 <- normalize_to_vcf_rep(
    chr = "chr1",
    pos = 502,
    ref = "GGTT",
    alt = "G_",
    fasta_fafile = fasta_loaded,
    one_based = TRUE,
    verbose = TRUE
  )
  expect_equal(out4$chr, "chr1")
  expect_equal(out4$pos, 502)
  expect_equal(out4$ref, "GGTT")
  expect_equal(out4$alt, "G")

  # --------------------------------------------------------------------------
  # Case 5: Complex deletion
  # --------------------------------------------------------------------------
  out5 <- normalize_to_vcf_rep(
    chr = "chr1",
    pos = 502,
    ref = "GGTT",
    alt = "A",
    fasta_fafile = fasta_loaded,
    one_based = TRUE,
    verbose = TRUE
  )
  expect_equal(out5$chr, "chr1")
  expect_equal(out5$pos, 502)
  expect_equal(out5$ref, "GGTT")
  expect_equal(out5$alt, "A")

  # --------------------------------------------------------------------------
  # Case 6: Complex insertion
  # --------------------------------------------------------------------------
  out6 <- normalize_to_vcf_rep(
    chr = "chr1",
    pos = 501,
    ref = "A",
    alt = "GAT",
    fasta_fafile = fasta_loaded,
    one_based = TRUE,
    verbose = TRUE
  )
  expect_equal(out6$chr, "chr1")
  expect_equal(out6$pos, 501)
  expect_equal(out6$ref, "A")
  expect_equal(out6$alt, "GAT")

  # --------------------------------------------------------------------------
  # Case 7: SNP, but position is given in 0-based coords
  # --------------------------------------------------------------------------
  out7 <- normalize_to_vcf_rep(
    chr = "chr1",
    pos = 502, # 0-based (matches 503 1-based)
    ref = "GT",
    alt = "AC",
    fasta_fafile = fasta_loaded,
    one_based = FALSE,
    verbose = TRUE
  )
  expect_equal(out7$chr, "chr1")
  expect_equal(out7$pos, 503)
  expect_equal(out7$ref, "GT")
  expect_equal(out7$alt, "AC")

  # --------------------------------------------------------------------------
  # Case 8: Check chromosome transformation
  # --------------------------------------------------------------------------
  out8 <- normalize_to_vcf_rep(
    chr = 1,
    pos = 503,
    ref = "G",
    alt = "A",
    fasta_fafile = fasta_loaded,
    one_based = TRUE,
    verbose = TRUE
  )
  expect_equal(out8$chr, "chr1")
  expect_equal(out8$pos, 503)
  expect_equal(out8$ref, "G")
  expect_equal(out8$alt, "A")

  # --------------------------------------------------------------------------
  # Case 9: Check if seq ref != fasta
  # --------------------------------------------------------------------------
  expect_warning(
    out9 <- normalize_to_vcf_rep(
      chr = "chr1", # Ensure this matches the FASTA convention
      pos = 503,
      ref = "T", # Incorrect reference allele (Actual is G)
      alt = "A",
      fasta_fafile = fasta_loaded,
      one_based = TRUE,
      verbose = TRUE
    ),
    "Mismatch found"
  )
  expect_null(out9)

  # --------------------------------------------------------------------------
  # Case 10: Expect no message when verbose = FALSE
  # --------------------------------------------------------------------------
  expect_no_message(
    out <- normalize_to_vcf_rep(
      chr = "chr1",
      pos = 503,
      ref = "G",
      alt = "A",
      fasta_fafile = fasta_loaded,
      one_based = TRUE,
      verbose = FALSE
    )
  )
  expect_equal(out$chr, "chr1")
  expect_equal(out$pos, 503)
  expect_equal(out$ref, "G")
  expect_equal(out$alt, "A")

  cleanup_test_fasta(fasta_env)
})
