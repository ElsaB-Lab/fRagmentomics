test_that("search_for_indel_in_cigar works for deletions", {
  # example 1 : deletion is correctly detected
  expect_equal(search_for_indel_in_cigar(
    pos        = 9,
    ref        = "TAT",
    alt        = "T",
    read_stats = list(POS=1, CIGAR="9M2D10M"),
    type       = "D"
  )[[1]], TRUE)

  # example 2: Case where the read doesn't carry the deletion
  expect_equal(search_for_indel_in_cigar(
    pos        = 9,
    ref        = "TA",
    alt        = "T",
    read_stats = list(POS=1, CIGAR="15M1D4M"),
    type       = "D"
  )[[1]], FALSE)

  # example 3: Case where the deletion length does not match the expected length
  expect_equal(search_for_indel_in_cigar(
    pos        = 9,
    ref        = "TAT",
    alt        = "T",
    read_stats = list(POS=1, CIGAR="9M1D5M"),
    type       = "D"
  )[[1]], FALSE)

  # example 4: Case where the read doesn't cover the position of the deletion
  expect_equal(search_for_indel_in_cigar(
    pos        = 9,
    ref        = "TAT",
    alt        = "T",
    read_stats = list(POS=1, CIGAR="8M"),
    type       = "D"
  )[[1]], FALSE)

  # example 5: Case where the read covers the position at the last nucleotide (nucleotide before the deletion)
  expect_equal(search_for_indel_in_cigar(
    pos        = 9,
    ref        = "TAT",
    alt        = "T",
    read_stats = list(POS=1, CIGAR="9M"),
    type       = "D"
  )[[1]], FALSE)


  # example 6: Case where there are 2 deletions and the 2nd one is correct
  expect_equal(search_for_indel_in_cigar(
    pos        = 9,
    ref        = "TAT",
    alt        = "T",
    read_stats = list(POS=1, CIGAR="1M2D6M2D10M"),
    type       = "D"
  )[[1]], TRUE)

  # example 7: Case where the read covers partially the deletion at the end
  expect_equal(search_for_indel_in_cigar(
    pos        = 9,
    ref        = "TAT",
    alt        = "T",
    read_stats = list(POS=1, CIGAR="9M1D"),
    type       = "D"
  )[[1]], FALSE)


  # example 8: Case where the read covers partially the deletion at the start
  expect_equal(search_for_indel_in_cigar(
    pos        = 1,
    ref        = "TAT",
    alt        = "T",
    read_stats = list(POS=1, CIGAR="1D9M"),
    type       = "D"
  )[[1]], FALSE)


  # example 9: Case where the read covers a larger deletion at the start
  expect_equal(search_for_indel_in_cigar(
    pos        = 1,
    ref        = "TAT",
    alt        = "T",
    read_stats = list(POS=1, CIGAR="3D9M"),
    type       = "D"
  )[[1]], FALSE)


  # example 10: Case where the read covers a larger deletion after the start
  expect_equal(search_for_indel_in_cigar(
    pos        = 1,
    ref        = "TAT",
    alt        = "T",
    read_stats = list(POS=1, CIGAR="1M3D9M"),
    type       = "D"
  )[[1]], FALSE)
})


test_that("search_for_indel_in_cigar works for insertions", {
  # example 1: Case where the insertion is correctly detected
  expect_equal(search_for_indel_in_cigar(
    pos        = 15,
    ref        = "A",
    alt        = "ATT",
    read_stats = list(POS=10, CIGAR="3S6M2I1M", SEQ="AAAAAAAAATTA"),
    type       = "I"
  )[[1]], TRUE)


  # example 2: Case where anoter insertion is present
  expect_equal(search_for_indel_in_cigar(
    pos        = 15,
    ref        = "A",
    alt        = "ATT",
    read_stats = list(POS=10, CIGAR="3S6M2I1M", SEQ="AAAAAAAAATCA"),
    type       = "I"
  )[[1]], FALSE)

  # example 3: Case where no insertion is present in the CIGAR
  expect_equal(search_for_indel_in_cigar(
    pos        = 15,
    ref        = "A",
    alt        = "ATT",
    read_stats = list(POS=10, CIGAR="3S7M", SEQ="AAAAAAAAAATT"),
    type       = "I"
  )[[1]], FALSE)

  # example 4: Case where the insertion length in CIGAR does not match alt
  expect_equal(search_for_indel_in_cigar(
    pos        = 15,
    ref        = "A",
    alt        = "ATT",
    read_stats = list(POS=10, CIGAR="3S6M1I2M", SEQ="AAAAAAAAAATT"),
    type       = "I"
  )[[1]], FALSE)

  # example 5: Case where the position does not correspond to the insertion area
  expect_equal(search_for_indel_in_cigar(
    pos        = 15,
    ref        = "A",
    alt        = "ATT",
    read_stats = list(POS=9, CIGAR="3S6M2I1M", SEQ="AAAAAAAAAATT"),
    type       = "I"
  )[[1]], FALSE)

  # example 6: Case with several insertions and the 2nd one good
  expect_equal(search_for_indel_in_cigar(
    pos        = 15,
    ref        = "A",
    alt        = "ATT",
    read_stats = list(POS=10, CIGAR="2M2I4M2I1M", SEQ="AAAAAAAATTA"),
    type       = "I"
  )[[1]], TRUE)

  # example 7: Case where the read covers completely the insertion at the start
  expect_equal(search_for_indel_in_cigar(
    pos        = 15,
    ref        = "A",
    alt        = "ATT",
    read_stats = list(POS=16, CIGAR="2I8M", SEQ="TTAAAAAAAA"),
    type       = "I"
  )[[1]], FALSE)

  # example 8: Case where another larger insertion is present at the start
  expect_equal(search_for_indel_in_cigar(
    pos        = 15,
    ref        = "A",
    alt        = "ATT",
    read_stats = list(POS=16, CIGAR="3I8M", SEQ="TTTAAAAAAAA"),
    type       = "I"
  )[[1]], FALSE)


  # example 9: Case where the read covers completely the insertion at the end
  expect_equal(search_for_indel_in_cigar(
    pos        = 3,
    ref        = "G",
    alt        = "GTT",
    read_stats = list(POS=1, CIGAR="3M2I", SEQ="AAGTT"),
    type       = "I"
  )[[1]], TRUE)

  # example 10: Case where the read covers completely the insertion but the nucleotide before the insertion is mutated
  expect_equal(search_for_indel_in_cigar(
    pos        = 3,
    ref        = "G",
    alt        = "GTT",
    read_stats = list(POS=1, CIGAR="3M2I2M", SEQ="AACTTCG"),
    type       = "I"
  )[[1]], FALSE)
})
