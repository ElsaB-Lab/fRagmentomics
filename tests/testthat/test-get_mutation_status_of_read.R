test_that("get_mutation_status_of_read works", {
  fasta <- system.file("testdata/fasta/hg19/", "hg19_chr1_27433000_27434000.fa", package = "fRagmentomics")
  fasta_fafile <- Rsamtools::FaFile(fasta)
  open(fasta_fafile)

  # INSERTIONS =============================================================================================

  # REF: AGTCCC
  # MUT: A > AGT
  # READ: AGTCC
  mstat <- get_mutation_status_of_read(
    chr                    = "chr1",
    pos                    = 500,
    ref                    = "A",
    alt                    = "AGT",
    read_index_at_pos      = 1,
    read_stats             = list(SEQ = "AGTCC", CIGAR = "5M", POS = 500),
    fasta_fafile           = fasta_fafile,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "WT")

  # REF: AGTCCC
  # MUT: A > AGT
  # READ: AGTGTC
  mstat <- get_mutation_status_of_read(
    chr                    = "chr1",
    pos                    = 500,
    ref                    = "A",
    alt                    = "AGT",
    read_index_at_pos      = 1,
    read_stats             = list(SEQ = "AGTGTC", CIGAR = "1M2I3M", POS = 500),
    fasta_fafile           = fasta_fafile,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "MUT")


  # REF: AGTCCC
  # MUT: A > AGT
  # READ: AGTGT
  mstat <- get_mutation_status_of_read(
    chr                    = "chr1",
    pos                    = 500,
    ref                    = "A",
    alt                    = "AGT",
    read_index_at_pos      = 1,
    read_stats             = list(SEQ = "AGTGT", CIGAR = "1M2I2M", POS = 500),
    fasta_fafile           = fasta_fafile,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "MUT by CIGAR but AMB")

  # Not the good position in the CIGAR
  mstat <- get_mutation_status_of_read(
    chr                    = "chr1",
    pos                    = 500,
    ref                    = "A",
    alt                    = "AGT",
    read_index_at_pos      = 1,
    read_stats             = list(SEQ = "AGTGT", CIGAR = "3M2I", POS = 500),
    fasta_fafile           = fasta_fafile,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "MUT but not in CIGAR and AMB")

  # REF: AGTCCC
  # MUT: A > AGT
  # READ: AGTG
  # CIGAR 3M1S
  mstat <- get_mutation_status_of_read(
    chr                    = "chr1",
    pos                    = 500,
    ref                    = "A",
    alt                    = "AGT",
    read_index_at_pos      = 1,
    read_stats             = list(SEQ = "AGTG", CIGAR = "3M1S", POS = 500),
    fasta_fafile           = fasta_fafile,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "MUT but not in CIGAR and AMB")

  # REF: AGTCCC
  # MUT: A > AGT
  # READ: AGT
  mstat <- get_mutation_status_of_read(
    chr                    = "chr1",
    pos                    = 500,
    ref                    = "A",
    alt                    = "AGT",
    read_index_at_pos      = 1,
    read_stats             = list(SEQ = "AGT", CIGAR = "3M", POS = 500),
    fasta_fafile           = fasta_fafile,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "AMB")

  # In CIGAR
  mstat <- get_mutation_status_of_read(
    chr                    = "chr1",
    pos                    = 500,
    ref                    = "A",
    alt                    = "AGT",
    read_index_at_pos      = 1,
    read_stats             = list(SEQ = "AGT", CIGAR = "1M2I", POS = 500),
    fasta_fafile           = fasta_fafile,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "MUT by CIGAR but AMB")

  # REF: AGTCCC
  # MUT: A > AGT
  # READ: AGTGAC
  mstat <- get_mutation_status_of_read(
    chr                    = "chr1",
    pos                    = 500,
    ref                    = "A",
    alt                    = "AGT",
    read_index_at_pos      = 1,
    read_stats             = list(SEQ = "AGTGAC", CIGAR = "3M2I1M", POS = 500),
    fasta_fafile           = fasta_fafile,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "OTH")

  # REF: GTCCCTCCTG
  # MUT: T > TC
  # READ: TCCCTC
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 502,
    ref = "T",
    alt = "TC",
    read_index_at_pos = 1,
    read_stats = list(SEQ = "TCCCTC", CIGAR = "6M", POS = 502),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "WT")

  # REF:  GTCCCTCCTG (501-510)
  # MUT: T > TC
  # READ: GTCCCCT
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 502,
    ref = "T",
    alt = "TC",
    read_index_at_pos = 2,
    read_stats = list(SEQ = "GTCCCCT", CIGAR = "2M1I4M", POS = 501),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT")

  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 502,
    ref = "T",
    alt = "TC",
    read_index_at_pos = 2,
    read_stats = list(SEQ = "GTCCCCT", CIGAR = "3M1I3M", POS = 501),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT but not in CIGAR")

  # REF:  GTCCCTCCTG (501-510)
  # MUT: T > TC
  # READ: GTCCCC
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 502,
    ref = "T",
    alt = "TC",
    read_index_at_pos = 2,
    read_stats = list(SEQ = "GTCCCC", CIGAR = "2M1I3M", POS = 501),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT by CIGAR but AMB")

  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 502,
    ref = "T",
    alt = "TC",
    read_index_at_pos = 2,
    read_stats = list(SEQ = "GTCCCC", CIGAR = "3M1I2M", POS = 501),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT but not in CIGAR and AMB")


  # DELETION COMPLETE COMPARISON =============================================================================================

  # REF TGGGAAGTCCCT (495-506)
  # MUT: GT > G (501)
  # READ: AAG
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 501,
    ref = "GT",
    alt = "G",
    read_index_at_pos = 3,
    read_stats = list(SEQ = "AAG", CIGAR = "3M", POS = 499),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "AMB")

  # REF TGGGAAGTCCCT (495-506)
  # MUT: GT > G (501)
  # READ: AAGCCC
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 501,
    ref = "GT",
    alt = "G",
    read_index_at_pos = 3,
    read_stats = list(SEQ = "AAGCCC", CIGAR = "3M1D3M", POS = 499),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT")


  # REF TGGGAAGTCCCT (495-506)
  # MUT: GT > G (501)
  # READ: AAGCC
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 501,
    ref = "GT",
    alt = "G",
    read_index_at_pos = 3,
    read_stats = list(SEQ = "AAGCC", CIGAR = "3M2D2M", POS = 499),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH by CIGAR but potentially MUT")

  # REF TGGGAAGTCCCT (495-506)
  # REF: GT / GC
  # MUT: GT > G (501)
  # READ: AAGTCC
  # Complex case with deletion and SNV
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 501,
    ref = "GT",
    alt = "G",
    read_index_at_pos = 3,
    read_stats = list(SEQ = "AAGTCC", CIGAR = "3M2D3M", POS = 499),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH by CIGAR but potentially WT")

  # REF TGGGAAGTCCCT (495-506)
  # MUT: AG > A (500)
  # READ: AAGTC
  # Complex case with deletion and SNV
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 500,
    ref = "AG",
    alt = "A",
    read_index_at_pos = 2,
    read_stats = list(SEQ = "AAGTC", CIGAR = "2M1D3M", POS = 499),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT by CIGAR but potentially WT")

  # REF TGGGAAGTCCCT (495-506)
  # MUT: AG > A (500)
  # READ: AATCC
  # Mutation position shifted
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 500,
    ref = "AG",
    alt = "A",
    read_index_at_pos = 2,
    read_stats = list(SEQ = "AATCC", CIGAR = "3M1D2M", POS = 499),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT but not in CIGAR")

  # REF TGGGAAGTCCCT (495-506)
  # REF: AG / MUT: AT
  # MUT: AG > A (500)
  # READ: AACCCT
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 500,
    ref = "AG",
    alt = "A",
    read_index_at_pos = 2,
    read_stats = list(SEQ = "AACCCT", CIGAR = "2M2D4M", POS = 499),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")

  # Deletion shifted to the right
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 500,
    ref = "AG",
    alt = "G",
    read_index_at_pos = 2,
    read_stats = list(SEQ = "AACCCT", CIGAR = "3M2D3M", POS = 499),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")

  # REF TGGGAAGTCCCT (495-506)
  # MUT: AG > A (500)
  # READ: AAACCT
  # SNV after the deletion
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 500,
    ref = "AG",
    alt = "G",
    read_index_at_pos = 2,
    read_stats = list(SEQ = "AAACCT", CIGAR = "2M1D4M", POS = 499),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT by CIGAR but potentially OTH")

  # REF TGGGAAGTCCCT (495-506)
  # MUT: AG > A (500)
  # READ: AAGGTCCCT
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 500,
    ref = "AG",
    alt = "G",
    read_index_at_pos = 2,
    read_stats = list(SEQ = "AAGTCCCT", CIGAR = "8M", POS = 499),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "WT")

  # DELETION COMPLETE COMPARISON =============================================================================================

  # REF TGGGAAGTCCCT (495-506)
  # REF:TGGG / MUT:TGGA
  # MUT: TG > T (495)
  # READ: TGG
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 495,
    ref = "TG",
    alt = "T",
    read_index_at_pos = 1,
    read_stats = list(SEQ = "TGG", CIGAR = "3M", POS = 495),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "AMB")

  # REF TGGGAAGTCCCT (495-506)
  # REF:TGGG / MUT:TGGA
  # MUT: TG > T (495)
  # READ: TGC
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 495,
    ref = "TG",
    alt = "T",
    read_index_at_pos = 1,
    read_stats = list(SEQ = "TGC", CIGAR = "3M", POS = 495),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")

  # REF TGGGAAGTCCCT (495-506)
  # REF:TGGG / MUT:TGGA
  # MUT: TG > T (495)
  # READ: TGG
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 495,
    ref = "TG",
    alt = "T",
    read_index_at_pos = 1,
    read_stats = list(SEQ = "TGG", CIGAR = "1M1D2M", POS = 495),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT by CIGAR but AMB")

  # REF TGGGAAGTCCCT (495-506)
  # REF:TGGG / MUT:TGGA
  # MUT: TG > T (495)
  # READ: TGC
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 495,
    ref = "TG",
    alt = "T",
    read_index_at_pos = 1,
    read_stats = list(SEQ = "TGC", CIGAR = "1M1D2M", POS = 495),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT by CIGAR but potentially OTH")

  # REF TGGGAAGTCCCT (495-506)
  # REF:TGGG / MUT:TGGA
  # MUT: TG > T (495)
  # READ: TGC
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 495,
    ref = "TG",
    alt = "T",
    read_index_at_pos = 1,
    read_stats = list(SEQ = "TGG", CIGAR = "1M2D2M", POS = 495),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH by CIGAR but AMB")

  # REF TGGGAAGTCCCT (495-506)
  # REF:TGGG / MUT:TGGA
  # MUT: TG > T (495)
  # READ: TGC
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 495,
    ref = "TG",
    alt = "T",
    read_index_at_pos = 1,
    read_stats = list(SEQ = "TGC", CIGAR = "1M2D2M", POS = 495),
    fasta_fafile = fasta_fafile,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")

  # COMPLEX CASES ======================================================================================================
  sequences <- c(chr1 = "ATCGAGGGGTCCAACCAAGGA")
  fasta_env <- setup_test_fasta(sequences)

  # REF:  ATCGAGGGGT
  # MUT: A > AGG
  # READ: ATCGAGGGGG
  mstat <- get_mutation_status_of_read(
    chr                    = "chr1",
    pos                    = 5,
    ref                    = "A",
    alt                    = "AGG",
    read_index_at_pos      = 5,
    read_stats             = list(SEQ = "ATCGAGGGGT", CIGAR = "11M", POS = 1),
    fasta_fafile           = fasta_env$fa_obj,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "WT")

  # REF:  ATCGAGGGGT
  # MUT: A > AGG
  # READ: ATCGAGGGGC
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGGC", CIGAR = "11M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")


  # REF:  ATCGAGGGGT
  # MUT: A > AGG
  # READ: ATCGAGGGGGG
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGGG", CIGAR = "5M1I4M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH by CIGAR but potentially MUT")

  # REF:  ATCGAGGGGT
  # MUT: A > AGG
  # READ: ATCGAGGGGTC
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGGTC", CIGAR = "5M1I5M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH by CIGAR but potentially WT")

  # REF:  ATCGAGGGGT
  # MUT: A > AGG
  # READ: ATCGAAGGGGG
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAAGGGGG", CIGAR = "5M1I5M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")

  # REF:  ATCGAGGGGT
  # MUT: A > AGG
  # READ: ATCGAGGGG
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGG", CIGAR = "5M1I4M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH by CIGAR but AMB")

  # REF:  ATCGAGGGGT
  # MUT: A > AGG
  # READ: ATCGAGGGT
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGT", CIGAR = "5M2I2M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT by CIGAR but potentially OTH")

  # REF:  ATCGAGGGGT
  # MUT: A > AGG
  # READ: ATCGAGGGG
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGG", CIGAR = "5M2I2M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT by CIGAR but AMB")

  # REF:  ATCGAGGGGT
  # MUT: A > AGG
  # READ: ATCGAGAGGGGT
  # Mutation at the end of the repeted sequence
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGGGAT", CIGAR = "5M2I5M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT by CIGAR but potentially OTH")

  # REF:  ATCGAGGGGT
  # MUT: A > AGG
  # READ: ATCGAGAGGGGT
  # Not the good inserted
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGAGGGGT", CIGAR = "5M2I5M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")

  # REF:  ATCGAGGGGT
  # MUT: A > AGG
  # READ: ATCGAGAGGGGT
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGACGGGGAT", CIGAR = "5M2I5M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")

  # REF:  ATCGAGGGGT
  # MUT: A > AGG
  # READ: ATCGAGAGGGGT
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGACGGGGGT", CIGAR = "5M2I5M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")

  # REF:  ATCGAGGGGT
  # MUT: A > AGG
  # READ: ATCGAGGGGGGGGT
  # insert of 4G instead of 2 at anchor position
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGGGGGGT", CIGAR = "5M4I3M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")

  # insert of 4G instead of 2 not at anchor position
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGGGGGGT", CIGAR = "7M4I1M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")

  cleanup_test_fasta(fasta_env)

  # complex case encountered in real-worl data =========================================================================

  sequences <- c(chr4 = "ACAGCACTATCTGAAACCAGGATGGATTGAATTGAAGGCCAAAGAGAGAGAAGAGATTTAGATGGATTTTAGAGTTCAAATGATATAG")
  fasta_env <- setup_test_fasta(sequences)

  # case 1: ATTGA from pos 26
  # ref_40 = 'ACAGCACTATCTGAAACCAGGATGGATTGAATTGAAGGCC'
  # seq_40 = 'ACAGCACTATCTGAAACCAGGATGG     ATTGAAGGCCCCTCG'
  # fake_ref= GATTGAAT / fake_alt= GATTGAAG
  mstat <- get_mutation_status_of_read(
    chr = "chr4",
    pos = 25,
    ref = "GATTGA",
    alt = "G",
    read_index_at_pos = 25,
    read_stats = list(SEQ = "ACAGCACTATCTGAAACCAGGATGGATTGAAGGCCCCTCG", CIGAR = "25M5D15M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT")

  # case 2: TTGAA from pos 32
  # ref_40 = 'ACAGCACTATCTGAAACCAGGATGGATTGAATTGAAGGCCAAAGAGAGAGAAGAGATTTAGATGGATTTTAGAGTTCAAATGATATAG'
  # seq_40 = 'ACAGCACTATCTGAAACCAGGATGGATTGAA     GGCCAAAGA'
  mstat <- get_mutation_status_of_read(
    chr = "chr4",
    pos = 31,
    ref = "ATTGAA",
    alt = "A",
    read_index_at_pos = 31,
    read_stats = list(SEQ = "ACAGCACTATCTGAAACCAGGATGGATTGAAGGCCAAAGA", CIGAR = "25M5D15M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT but not in CIGAR")

  # Case 3: Cigar with 2 deletions inside with the same size
  # ref_40 = 'ACAGCACTATCTGAAACCAGGATGGATTGAATTGAAGGCCAAAGAGAGAGAAGAGATTTAGATGGATTTTAGAGTTCAAATGATATAG'
  # ref_35 = 'A     CTATCTGAAACCAGGATGGATTGA     AGGCCAAAGA'
  mstat <- get_mutation_status_of_read(
    chr = "chr4",
    pos = 30,
    ref = "AATTGA",
    alt = "A",
    read_index_at_pos = 25,
    read_stats = list(SEQ = "ACTATCTGAAACCAGGATGGATTGAAGGCCAAAGA", CIGAR = "1M5D24M5D10M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT")

  # Case 4: OTH because softclipping is considered as real sequence
  # ref_40 = 'ACAGCACTATCTGAAACCAGGATGGATTGAATTGAAGGCC'
  # seq_31 = 'ACAGCACTATCTGAAACCAGGSSSSSSSSSS'
  # fake_ref = GA  / fake_mut = GG
  mstat <- get_mutation_status_of_read(
    chr = "chr4",
    pos = 21,
    ref = "GAT",
    alt = "G",
    read_index_at_pos = 21,
    read_stats = list(SEQ = "ACAGCACTATCTGAAACCAGGTTTTTTTTTT", CIGAR = "21M10S", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")

  # Case 5: MNV when a part is in soft clipping
  # ref_40 = 'ACAGCACTATCTGAAACCAGGATGGATTGAATTGAAGGCC'
  # seq_31 = 'ACAGCACTATCTGAAACCACCCTTTTTTTTT'
  mstat <- get_mutation_status_of_read(
    chr = "chr4",
    pos = 20,
    ref = "GGA",
    alt = "CCC",
    read_index_at_pos = 20,
    read_stats = list(SEQ = "ACAGCACTATCTGAAACCACCCTTTTTTTTT", CIGAR = "21M10S", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT")

  # Case 6: MNV when a part is absent
  # ref_40 = 'ACAGCACTATCTGAAACCAGGATGGATTGAATTGAAGGCC'
  # seq_21 = 'ACAGCACTATCTGAAACCACC'
  mstat <- get_mutation_status_of_read(
    chr = "chr4",
    pos = 20,
    ref = "GGA",
    alt = "CCC",
    read_index_at_pos = 20,
    read_stats = list(SEQ = "ACAGCACTATCTGAAACCACC", CIGAR = "21M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "AMB")

  # Case 8: MNV when a part is in soft clipping matching the ref
  # ref_40 = 'ACAGCACTATCTGAAACCAGGATGGATTGAATTGAAGGCC'
  # seq_21 = 'ACAGCACTATCTGAAACCACCCT'
  mstat <- get_mutation_status_of_read(
    chr = "chr4",
    pos = 20,
    ref = "GGA",
    alt = "CCC",
    read_index_at_pos = 20,
    read_stats = list(SEQ = "ACAGCACTATCTGAAACCACCCT", CIGAR = "21M2S", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT")

  # Case 9: MNV when a part is in soft clipping not matching the ref
  # ref_40 = 'ACAGCACTATCTGAAACCAGGATGGATTGAATTGAAGGCC'
  # seq_21 = 'ACAGCACTATCTGAAACCACCCAAA'
  mstat <- get_mutation_status_of_read(
    chr = "chr4",
    pos = 20,
    ref = "GGA",
    alt = "CCC",
    read_index_at_pos = 20,
    read_stats = list(SEQ = "ACAGCACTATCTGAAACCACCCAAA", CIGAR = "21M4S", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")

  # Case 10: Other Mut in the condition of AMB
  # ref_40 = 'ACAGCACTATCTGAAACCAGGATGGATTGAATTGAAGGCC'
  # seq_20 = 'ACAGCACTATCTGAAACCAC'
  mstat <- get_mutation_status_of_read(
    chr = "chr4",
    pos = 19,
    ref = "AG",
    alt = "A",
    read_index_at_pos = 19,
    read_stats = list(SEQ = "ACAGCACTATCTGAAACCAC", CIGAR = "20M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")

  cleanup_test_fasta(fasta_env)
})
