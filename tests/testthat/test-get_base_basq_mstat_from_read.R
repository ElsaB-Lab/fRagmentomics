test_that("get_base_basq_mstat_from_read works", {
  sequences <- c(chr1 = "ATCGAGGGGTCCAACCAAGGA")
  fasta_env <- setup_test_fasta(sequences)

  # Case where the position of the mutation is not covered
  mstat <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 50,
    ref = "A",
    alt = "AGG",
    read_stats = list(SEQ = "ATCGAGGGGG", QUAL = "##########", CIGAR = "10M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_mode = FALSE
  )
  expect_equal(mstat, list(base = NA, basq = NA, mstat = NA))

  # Case where we cannot get back all the necessary sequence
  mstat <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_stats = list(SEQ = "ATCGAG", QUAL = "####!!", CIGAR = "6M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_mode = FALSE
  )
  expect_equal(mstat, list(base = "AG*", basq = "!!*", mstat = "AMB"))

  # Normal SNV case
  mstat <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 7,
    ref = "A",
    alt = "T",
    read_stats = list(SEQ = "ATCGTGG", QUAL = "####!##", CIGAR = "1M2D6M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_mode = FALSE
  )
  expect_equal(mstat, list(base = "T", basq = "!", mstat = "MUT"))

  # Normal indel case
  mstat <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_stats = list(SEQ = "ATCGAGGGGT", QUAL = "####!!!###", CIGAR = "5M2I3M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_mode = FALSE
  )
  expect_equal(mstat, list(base = "AGG", basq = "!!!", mstat = "MUT"))

  cleanup_test_fasta(fasta_env)
})
