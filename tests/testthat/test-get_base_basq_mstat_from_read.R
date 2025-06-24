test_that("get_base_basq_mstat_from_read works", {
  sequences <- c(chr1 = "ATCGAGGGGTCCAACCAAGGA")
  fasta_env <- setup_test_fasta(sequences)

  # Case where the position of the mutation is not covered
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 50,
    ref = "A",
    alt = "AGG",
    read_stats = list(SEQ = "ATCGAGGGGG", QUAL = "##########", CIGAR = "10M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE
  )
  expect_equal(read_info, list(base = NA, basq = NA, mstat = NA))

  # Case where we cannot get back all the necessary sequence
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_stats = list(SEQ = "ATCGAG", QUAL = "####!!", CIGAR = "6M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE
  )
  expect_equal(read_info, list(base = "AG*", basq = "!!*", mstat = "AMB"))

  # Normal SNV case
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 7,
    ref = "A",
    alt = "T",
    read_stats = list(SEQ = "ATCGTGG", QUAL = "####!##", CIGAR = "1M2D6M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE
  )
  expect_equal(read_info, list(base = "T", basq = "!", mstat = "MUT"))

  # Normal INS case
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_stats = list(SEQ = "ATCGAGGGGT", QUAL = "####!!!###", CIGAR = "5M2I3M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE
  )
  expect_equal(read_info, list(base = "AGG", basq = "!!!", mstat = "MUT"))


  # SNV in soft-clipping
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 6,
    ref = "G",
    alt = "T",
    read_stats = list(SEQ = "ATCGATTTTA", QUAL = "#A!$Q###A!", CIGAR = "5M5S", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE
  )
  expect_equal(read_info, list(base = NA, basq = NA, mstat = NA))

  # Normal MNV case
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 6,
    ref = "GG",
    alt = "CC",
    read_stats = list(SEQ = "ATCGACCGGT", QUAL = "#A!$Q###A!", CIGAR = "10M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE
  )
  expect_equal(read_info, list(base = "CC", basq = "##", mstat = "MUT"))

  # Other MNV case
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 6,
    ref = "GG",
    alt = "CC",
    read_stats = list(SEQ = "ATCGACAGGT", QUAL = "#A!$Q###A!", CIGAR = "10M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE
  )
  expect_equal(read_info, list(base = "CA", basq = "##", mstat = "OTH"))


  # SNV flanked by another SNV
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 6,
    ref = "G",
    alt = "C",
    read_stats = list(SEQ = "ATCGACAGGT", QUAL = "#A!$Q###A!", CIGAR = "10M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE
  )
  expect_equal(read_info, list(base = "C", basq = "#", mstat = "MUT but potentially larger MNV"))

  # MNV flanked by another SNV
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 6,
    ref = "GG",
    alt = "CC",
    read_stats = list(SEQ = "ATCGACCCGT", QUAL = "#A!$Q###A!", CIGAR = "10M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE
  )
  expect_equal(read_info, list(base = "CC", basq = "##", mstat = "MUT but potentially larger MNV"))


  # SNV in DEL case
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 10,
    ref = "T",
    alt = "G",
    read_stats = list(SEQ = "ATCGAGGGGACCAA", QUAL = "#A!$Q#UAOQUELA", CIGAR = "9M4D5M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE
  )
  expect_equal(read_info, list(base = "-", basq = "", mstat = "OTH (DEL)"))

  # SNV case in the border of an indel
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 4,
    ref = "G",
    alt = "T",
    read_stats = list(SEQ = "ATAGGGG", QUAL = "#!#####", CIGAR = "1M2D6M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE
  )
  expect_equal(read_info, list(base = "T", basq = "!", mstat = "MUT but potentially larger MNV"))

  # SNV case at the beginning of the read
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 1,
    ref = "A",
    alt = "T",
    read_stats = list(SEQ = "TTCGTGG", QUAL = "!######", CIGAR = "1M2D6M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE
  )
  expect_equal(read_info, list(base = "T", basq = "!", mstat = "MUT"))

  # SNV case in the border of an indel (ins) and last nucleotide
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "G",
    read_stats = list(SEQ = "ATCGTTG", QUAL = "######!", CIGAR = "4M2I1M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE
  )
  expect_equal(read_info, list(base = "G", basq = "!", mstat = "MUT but potentially larger MNV"))

  # "ATCGAGGGGTCCAACCAAGGA"

  cleanup_test_fasta(fasta_env)
})
