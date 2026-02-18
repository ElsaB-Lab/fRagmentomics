# Define the toy sequence and setup the FASTA environment for all tests
sequences <- c(chr4 = "ACAGCACTATCTGAAACCAGGATGGATTGAATTGAAGGCCAAAGAGAGAGAAGAGATTTAGATGGATTTTAGAGTTCAAATGATATAG")
fasta_env <- setup_test_fasta(sequences)
fasta_fafile <- fasta_env$fa_obj

# =========================================================================
## Tests for get_seq_from_fasta()
# =========================================================================

test_that("get_seq_from_fasta works correctly with a pre-loaded fasta_seq object", {
  fasta_seq <- list(
    chr = "chr4",
    start = 10,
    end = 30,
    seq = "TCTGAAACCAGGATGGATTGA"
  )

  # Test 1: Successful extraction from the middle of the pre-loaded sequence
  expect_equal(
    get_seq_from_fasta(chr = "chr4", start = 15, end = 20, fasta_seq = fasta_seq),
    "AACCAG"
  )

  # Test 2: Error when requested chromosome does not match
  expect_error(
    get_seq_from_fasta(chr = "chr5", start = 15, end = 20, fasta_seq = fasta_seq),
    "Requested chromosome 'chr5' does not match the available chromosome 'chr4'."
  )

  # Test 3: Error when requested range is outside the available sequence boundaries
  expect_error(
    get_seq_from_fasta(chr = "chr4", start = 5, end = 15, fasta_seq = fasta_seq),
    "Requested sequence range 5:15 does not fit into the available reference sequence range 10:30."
  )
  expect_error(
    get_seq_from_fasta(chr = "chr4", start = 25, end = 35, fasta_seq = fasta_seq),
    "Requested sequence range 25:35 does not fit into the available reference sequence range 10:30."
  )
})

test_that("get_seq_from_fasta works correctly with a FaFile object", {
  # Test 1: Fetch a single nucleotide
  expect_equal(
    get_seq_from_fasta(chr = "chr4", start = 1, end = 1, fasta_fafile = fasta_fafile),
    "A"
  )

  # Test 2: Fetch a sequence range
  expect_equal(
    get_seq_from_fasta(chr = "chr4", start = 1, end = 10, fasta_fafile = fasta_fafile),
    "ACAGCACTAT"
  )
})


# =========================================================================
## Tests for harmonize_chr_to_fasta()
# =========================================================================

test_that("harmonize_chr_to_fasta correctly adds or removes 'chr' prefix", {
  # --- Test Case 1: FASTA uses 'chr' prefix (using the global fasta_fafile) ---
  # Adds 'chr' prefix when missing
  expect_equal(harmonize_chr_to_fasta("4", fasta_fafile), "chr4")
  # Does nothing if 'chr' prefix is already present
  expect_equal(harmonize_chr_to_fasta("chr4", fasta_fafile), "chr4")

  # --- Test Case 2: FASTA does not use 'chr' prefix ---
  # Create a temporary FASTA file without 'chr' prefixes for this test
  sequences_no_chr <- c("4" = "ACAGCACTAT")
  fasta_env_no_chr <- setup_test_fasta(sequences_no_chr)
  fasta_fafile_no_chr <- fasta_env_no_chr$fa_obj

  # Removes 'chr' prefix when present
  expect_equal(harmonize_chr_to_fasta("chr4", fasta_fafile_no_chr), "4")
  # Does nothing if 'chr' prefix is already absent
  expect_equal(harmonize_chr_to_fasta("4", fasta_fafile_no_chr), "4")

  # Clean up temporary file for the no-chr case
  cleanup_test_fasta(fasta_env_no_chr)
})


# =========================================================================
## Tests for check_if_ref_matches_fasta()
# =========================================================================

test_that("check_if_ref_matches_fasta correctly validates reference alleles", {
  # Test 1: Correct single nucleotide reference returns TRUE
  expect_true(check_if_ref_matches_fasta(chr = "chr4", pos = 1, ref = "A", fasta_fafile = fasta_fafile))

  # Test 2: Correct multi-nucleotide reference (insertion) returns TRUE
  expect_true(check_if_ref_matches_fasta(chr = "chr4", pos = 10, ref = "TCTGAA", fasta_fafile = fasta_fafile))

  # Test 3: Incorrect single nucleotide reference returns FALSE
  expect_false(check_if_ref_matches_fasta(chr = "chr4", pos = 1, ref = "C", fasta_fafile = fasta_fafile))

  # Test 4: Incorrect multi-nucleotide reference returns FALSE
  expect_false(check_if_ref_matches_fasta(chr = "chr4", pos = 10, ref = "AAAAAA", fasta_fafile = fasta_fafile))

  # Test 5: Chromosome not found in FASTA issues a warning and returns FALSE
  expect_warning(
    result <- check_if_ref_matches_fasta(chr = "chrX", pos = 1, ref = "A", fasta_fafile = fasta_fafile),
    "Chromosome 'chrX' not found in FASTA."
  )
  expect_false(result)
})

cleanup_test_fasta(fasta_env)
