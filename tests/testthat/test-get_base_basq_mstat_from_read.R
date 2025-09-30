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
  )
  expect_equal(read_info, list(base = "AGG", basq = "!!!", mstat = "MUT by CIGAR but potentially WT"))

  # WT INS case without repetition
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 1,
    ref = "A",
    alt = "AGG",
    read_stats = list(SEQ = "ATCGAGGGGT", QUAL = "abcdefghij", CIGAR = "5M2I3M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
  )
  expect_equal(read_info, list(base = "ATC", basq = "abc", mstat = "WT"))

  # WT DEL case without repetition
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 1,
    ref = "ATC",
    alt = "A",
    read_stats = list(SEQ = "ATCGAGGGGT", QUAL = "*###!!!###", CIGAR = "5M2I3M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
  )
  expect_equal(read_info, list(base = "A", basq = "*", mstat = "WT"))

  # SNV in soft-clipping
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 6,
    ref = "G",
    alt = "T",
    read_stats = list(SEQ = "ATCGATTTTA", QUAL = "abcdefghij", CIGAR = "5M5S", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
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
  )
  expect_equal(read_info, list(base = "C", basq = "#", mstat = "MUT but potentially OTH"))

  # MNV flanked by another SNV
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 6,
    ref = "GG",
    alt = "CC",
    read_stats = list(SEQ = "ATCGACCCGT", QUAL = "#A!$Q###A!", CIGAR = "10M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
  )
  expect_equal(read_info, list(base = "CC", basq = "##", mstat = "MUT but potentially OTH"))

  # MNV case where a part of the MNV is in softclipped bases
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 6,
    ref = "GG",
    alt = "CC",
    read_stats = list(SEQ = "ATCGACCGGT", QUAL = "#A!$Q###A!", CIGAR = "6M4S", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
  )
  # TODO: expect_equal(read_info, list(base = "C*", basq = "#*", mstat = "AMB"))
  expect_equal(read_info, list(base = "CC", basq = "##", mstat = "MUT"))

  # Del case where the anchor position is before softclipping
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 5,
    ref = "AGGGG",
    alt = "A",
    read_stats = list(SEQ = "ATCGATCCAA", QUAL = "#A!$Q###A!", CIGAR = "5M5S", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
  )
  expect_equal(read_info, list(base = "A", basq = "Q", mstat = "MUT but not in CIGAR"))

  # Del case where the anchor position is in softclipping
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 5,
    ref = "AGGGG",
    alt = "A",
    read_stats = list(SEQ = "ATCGATCCAA", QUAL = "#A!$Q###A!", CIGAR = "4M6S", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
  )
  expect_equal(read_info, list(base = NA, basq = NA, mstat = NA))

  # Ins case where the anchor position is before softclipping
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 6,
    ref = "G",
    alt = "GTCA",
    read_stats = list(SEQ = "ATCGAGTCAG", QUAL = "#A!$Q###A!", CIGAR = "6M4S", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
  )
  expect_equal(read_info, list(base = "GTCA", basq = "###A", mstat = "MUT but not in CIGAR"))

  # Ins case where the anchor position is in softclipping
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 6,
    ref = "G",
    alt = "GTCA",
    read_stats = list(SEQ = "ATCGAGTCAG", QUAL = "#A!$Q###A!", CIGAR = "5M5S", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
  )
  expect_equal(read_info, list(base = NA, basq = NA, mstat = NA))

  # SNV in DEL case
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 10,
    ref = "T",
    alt = "G",
    read_stats = list(SEQ = "ATCGAGGGGACCAA", QUAL = "#A!$Q#UAOQUELA", CIGAR = "9M4D5M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
  )
  expect_equal(read_info, list(base = "-", basq = "", mstat = "OTH"))

  # SNV case in the border of an indel
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 4,
    ref = "G",
    alt = "T",
    read_stats = list(SEQ = "ATAGGGG", QUAL = "#!#####", CIGAR = "1M2D6M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
  )
  expect_equal(read_info, list(base = "T", basq = "!", mstat = "MUT but potentially OTH"))

  # SNV case at the beginning of the read
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 1,
    ref = "A",
    alt = "T",
    read_stats = list(SEQ = "TTCGTGG", QUAL = "!######", CIGAR = "1M2D6M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
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
  )
  expect_equal(read_info, list(base = "G", basq = "!", mstat = "MUT but potentially OTH"))

  cleanup_test_fasta(fasta_env)

  # new sequence
  sequences <- c(chr1 = "ACTGTGTGTGTGATATATGCGAGAGAGT")
  fasta_env <- setup_test_fasta(sequences)

  # REF = TGTGTGTGTGA
  # TGTGTGTGTGT
  # MUT = TGTGTGTGTGTGTGA
  # TGTGTGTGTGTGTGA
  # CIGAR reports insert of GT but we look for insertion of GTGT
  read_info <- get_base_basq_mstat_from_read(
    chr = "chr1",
    pos = 3,
    ref = "T",
    alt = "TGTGT",
    read_stats = list(SEQ = "ACTGTGTGTGTGTGTGA", QUAL = "abcdefghij#!####B", CIGAR = "3M2I12M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
  )
  expect_equal(read_info, list(base = "TGTGT", basq = "cdefg", mstat = "OTH by CIGAR but potentially MUT"))


  cleanup_test_fasta(fasta_env)
})
