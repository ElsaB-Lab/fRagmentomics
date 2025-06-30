test_that("get_mutation_status_of_read works", {
  fasta <- system.file("testdata/fasta/hg19/", "hg19_chr1_27433000_27434000.fa", package = "fRagmentomics")
  fasta_fafile <- Rsamtools::FaFile(fasta)
  open(fasta_fafile)

  # INSERTIONS NO REPEAT =============================================================================================

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
    cigar_free_indel_match = FALSE,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "WT")

  mstat <- get_mutation_status_of_read(
    chr                    = "chr1",
    pos                    = 500,
    ref                    = "A",
    alt                    = "AGT",
    read_index_at_pos      = 1,
    read_stats             = list(SEQ = "AGTCC", CIGAR = "5M", POS = 500),
    fasta_fafile           = fasta_fafile,
    cigar_free_indel_match = TRUE,
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
    cigar_free_indel_match = FALSE,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "MUT")


  mstat <- get_mutation_status_of_read(
    chr                    = "chr1",
    pos                    = 500,
    ref                    = "A",
    alt                    = "AGT",
    read_index_at_pos      = 1,
    read_stats             = list(SEQ = "AGTGTC", CIGAR = "1M2I3M", POS = 500),
    fasta_fafile           = fasta_fafile,
    cigar_free_indel_match = TRUE,
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
    cigar_free_indel_match = FALSE,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "MUT")

  mstat <- get_mutation_status_of_read(
    chr                    = "chr1",
    pos                    = 500,
    ref                    = "A",
    alt                    = "AGT",
    read_index_at_pos      = 1,
    read_stats             = list(SEQ = "AGTGT", CIGAR = "3M2I", POS = 500),
    fasta_fafile           = fasta_fafile,
    cigar_free_indel_match = FALSE,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "AMB")

  mstat <- get_mutation_status_of_read(
    chr                    = "chr1",
    pos                    = 500,
    ref                    = "A",
    alt                    = "AGT",
    read_index_at_pos      = 1,
    read_stats             = list(SEQ = "AGTGT", CIGAR = "3M2I", POS = 500),
    fasta_fafile           = fasta_fafile,
    cigar_free_indel_match = TRUE,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "AMB by cigar-free search and MUT not found in CIGAR")

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
    cigar_free_indel_match = FALSE,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "AMB")

  mstat <- get_mutation_status_of_read(
    chr                    = "chr1",
    pos                    = 500,
    ref                    = "A",
    alt                    = "AGT",
    read_index_at_pos      = 1,
    read_stats             = list(SEQ = "AGT", CIGAR = "1M2I", POS = 500),
    fasta_fafile           = fasta_fafile,
    cigar_free_indel_match = FALSE,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "MUT")

  mstat <- get_mutation_status_of_read(
    chr                    = "chr1",
    pos                    = 500,
    ref                    = "A",
    alt                    = "AGT",
    read_index_at_pos      = 1,
    read_stats             = list(SEQ = "AGT", CIGAR = "3M", POS = 500),
    fasta_fafile           = fasta_fafile,
    cigar_free_indel_match = TRUE,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "AMB by cigar-free search and MUT not found in CIGAR")

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
    cigar_free_indel_match = FALSE,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "WT")

  mstat <- get_mutation_status_of_read(
    chr                    = "chr1",
    pos                    = 500,
    ref                    = "A",
    alt                    = "AGT",
    read_index_at_pos      = 1,
    read_stats             = list(SEQ = "AGTGAC", CIGAR = "3M2I1M", POS = 500),
    fasta_fafile           = fasta_fafile,
    cigar_free_indel_match = TRUE,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "OTH")

  # INSERTIONS REPEAT ================================================================================================

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
    cigar_free_indel_match = TRUE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "WT")

  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 502,
    ref = "T",
    alt = "TC",
    read_index_at_pos = 1,
    read_stats = list(SEQ = "TCCCTC", CIGAR = "6M", POS = 502),
    fasta_fafile = fasta_fafile,
    cigar_free_indel_match = TRUE,
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
    cigar_free_indel_match = TRUE,
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
    cigar_free_indel_match = TRUE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT by cigar-free search and MUT not found in CIGAR")


  # DELETION NO REPEAT =============================================================================================

  # REF TGGGAAGTCCCT (495-506)
  # MUT: GT > T (501)
  # READ: AAGCCC
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 501,
    ref = "GT",
    alt = "T",
    read_index_at_pos = 3,
    read_stats = list(SEQ = "AAGCCC", CIGAR = "3M1D3M", POS = 499),
    fasta_fafile = fasta_fafile,
    cigar_free_indel_match = TRUE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT")

  # REF TGGGAAGTCCCT (495-506)
  # MUT: GT > T (501)
  # READ: AAGCC
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 501,
    ref = "GT",
    alt = "T",
    read_index_at_pos = 3,
    read_stats = list(SEQ = "AAGCC", CIGAR = "3M2D2M", POS = 499),
    fasta_fafile = fasta_fafile,
    cigar_free_indel_match = FALSE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")

  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 501,
    ref = "GT",
    alt = "T",
    read_index_at_pos = 3,
    read_stats = list(SEQ = "AAGCC", CIGAR = "3M2D2M", POS = 499),
    fasta_fafile = fasta_fafile,
    cigar_free_indel_match = TRUE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT by cigar-free search and other MUT found in CIGAR")


  # REF TGGGAAGTCCCT (495-506)
  # MUT: GT > T (501)
  # READ: AAGCCT
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 501,
    ref = "GT",
    alt = "T",
    read_index_at_pos = 3,
    read_stats = list(SEQ = "AAGCCT", CIGAR = "3M2D3M", POS = 499),
    fasta_fafile = fasta_fafile,
    cigar_free_indel_match = FALSE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")

  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 501,
    ref = "GT",
    alt = "T",
    read_index_at_pos = 3,
    read_stats = list(SEQ = "AAGCCT", CIGAR = "3M2D3M", POS = 499),
    fasta_fafile = fasta_fafile,
    cigar_free_indel_match = TRUE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT by cigar-free search and other MUT found in CIGAR")


  # REF TGGGAAGTCCCT (495-506)
  # MUT: AG > A (500)
  # READ: AATCC
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 500,
    ref = "AG",
    alt = "A",
    read_index_at_pos = 2,
    read_stats = list(SEQ = "AATCC", CIGAR = "2M1D3M", POS = 499),
    fasta_fafile = fasta_fafile,
    cigar_free_indel_match = FALSE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT")

  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 500,
    ref = "AG",
    alt = "A",
    read_index_at_pos = 2,
    read_stats = list(SEQ = "AATCC", CIGAR = "2M1D3M", POS = 499),
    fasta_fafile = fasta_fafile,
    cigar_free_indel_match = TRUE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT")


  # REF TGGGAAGTCCCT (495-506)
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
    cigar_free_indel_match = FALSE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")

  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 500,
    ref = "AG",
    alt = "G",
    read_index_at_pos = 2,
    read_stats = list(SEQ = "AACCCT", CIGAR = "2M2D4M", POS = 499),
    fasta_fafile = fasta_fafile,
    cigar_free_indel_match = TRUE,
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
    read_stats             = list(SEQ = "ATCGAGGGGG", CIGAR = "11M", POS = 1),
    fasta_fafile           = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE,
    n_match_base_before    = 1,
    n_match_base_after     = 1
  )
  expect_equal(mstat, "AMB")

  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGGG", CIGAR = "11M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = TRUE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "AMB by cigar-free search and MUT not found in CIGAR")


  # REF:  ATCGAGGGGT
  # MUT: A > AGG
  # READ: ATCGAGGGGGG
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGGGG", CIGAR = "12M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "AMB")

  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGGGG", CIGAR = "12M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = TRUE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "AMB by cigar-free search and MUT not found in CIGAR")

  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGGGG", CIGAR = "5M2I4M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT")

  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGGGG", CIGAR = "7M2I2M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "AMB")

  # REF:  ATCGAGGGGT
  # MUT: A > AGG
  # READ: ATCGAGGGGGGT
  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGGGGT", CIGAR = "5M2I5M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT")

  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGGGGT", CIGAR = "5M2I5M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = TRUE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT")

  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGGGGT", CIGAR = "7M2I3M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "WT")

  mstat <- get_mutation_status_of_read(
    chr = "chr1",
    pos = 5,
    ref = "A",
    alt = "AGG",
    read_index_at_pos = 5,
    read_stats = list(SEQ = "ATCGAGGGGGGT", CIGAR = "7M2I3M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = TRUE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT by cigar-free search and MUT not found in CIGAR")

  cleanup_test_fasta(fasta_env)

  sequences <- c(chr4 = "ACAGCACTATCTGAAACCAGGATGGATTGAATTGAAGGCCAAAGAGAGAGAAGAGATTTAGATGGATTTTAGAGTTCAAATGATATAG")
  fasta_env <- setup_test_fasta(sequences)

  # case 1: ATTGA from pos 26
  # ref_40 = 'ACAGCACTATCTGAAACCAGGATGGATTGAATTGAAGGCC'
  # seq_40 = 'ACAGCACTATCTGAAACCAGGATGG     ATTGAAGGCCCCTCG'
  mstat <- get_mutation_status_of_read(
    chr = "chr4",
    pos = 25,
    ref = "GATTGA",
    alt = "G",
    read_index_at_pos = 25,
    read_stats = list(SEQ = "ACAGCACTATCTGAAACCAGGATGGATTGAAGGCCCCTCG", CIGAR = "25M5D15M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT")

  mstat <- get_mutation_status_of_read(
    chr = "chr4",
    pos = 25,
    ref = "GATTGA",
    alt = "G",
    read_index_at_pos = 25,
    read_stats = list(SEQ = "ACAGCACTATCTGAAACCAGGATGGATTGAAGGCCCCTCG", CIGAR = "25M5D15M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = TRUE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT")

  # case 2: TTGAA from pos 32
  # ref_40 = 'ACAGCACTATCTGAAACCAGGATGGATTGAATTGAAGGCC'
  # seq_40 = 'ACAGCACTATCTGAAACCAGGATGGATTGAA     GGCCCCTCG'
  mstat <- get_mutation_status_of_read(
    chr = "chr4",
    pos = 31,
    ref = "ATTGAA",
    alt = "A",
    read_index_at_pos = 31,
    read_stats = list(SEQ = "ACAGCACTATCTGAAACCAGGATGGATTGAAGGCCCCTCG", CIGAR = "25M5D15M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "WT")

  mstat <- get_mutation_status_of_read(
    chr = "chr4",
    pos = 31,
    ref = "ATTGAA",
    alt = "A",
    read_index_at_pos = 31,
    read_stats = list(SEQ = "ACAGCACTATCTGAAACCAGGATGGATTGAAGGCCCCTCG", CIGAR = "25M5D15M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = TRUE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT by cigar-free search and MUT not found in CIGAR")


  # Case 3: Cigar with 2 deletions inside with the same size
  mstat <- get_mutation_status_of_read(
    chr = "chr4",
    pos = 30,
    ref = "ATTGAA",
    alt = "A",
    read_index_at_pos = 25,
    read_stats = list(SEQ = "ACAGCACTATCTGAAACCAGGATGGATTGAAGGCCCCTCG", CIGAR = "1M5D24M5D15M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "MUT")

  # Case 4: Ambiguous when the position is before softclipping
  mstat <- get_mutation_status_of_read(
    chr = "chr4",
    pos = 21,
    ref = "GAT",
    alt = "G",
    read_index_at_pos = 21,
    read_stats = list(SEQ = "ACAGCACTATCTGAAACCAGG", CIGAR = "21M10S", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "AMB")

  # Case 5: MNV when a part is in soft clipping
  mstat <- get_mutation_status_of_read(
    chr = "chr4",
    pos = 20,
    ref = "GGA",
    alt = "CCC",
    read_index_at_pos = 20,
    read_stats = list(SEQ = "ACAGCACTATCTGAAACCACC", CIGAR = "21M10S", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "AMB")

  # Case 6: Other Mut in the condition of AMB
  mstat <- get_mutation_status_of_read(
    chr = "chr4",
    pos = 19,
    ref = "AG",
    alt = "A",
    read_index_at_pos = 19,
    read_stats = list(SEQ = "ACAGCACTATCTGAAACCAC", CIGAR = "20M", POS = 1),
    fasta_fafile = fasta_env$fa_obj,
    cigar_free_indel_match = FALSE,
    n_match_base_before = 1,
    n_match_base_after = 1
  )
  expect_equal(mstat, "OTH")

  cleanup_test_fasta(fasta_env)
})
