test_that("get insertion", {
  # example 1: Case where the insertion is correctly detected

  res1 <- get_insertion(
    pos      = 15,
    alt      = "ATT",
    r_pos    = 10,
    r_cigar  = "3S6M2I1M",
    r_query  = "AAAAAAAAATTA",
    r_qual   = "########I++#"
  )

  expect_equal(res1$base, "+TT")
  expect_equal(res1$qual, "I")

  # example 2: Case where no insertion is present in the CIGAR
  res2 <- get_insertion(
    pos      = 15,
    alt      = "ATT",
    r_pos    = 10,
    r_cigar  = "3S7M",
    r_query  = "AAAAAAAAAATT",
    r_qual   = "##########++"
  )

  expect_equal(res2$base, "no_insertion_detected")
  expect_equal(res2$qual, "no_insertion_detected")


  # example 3: Case where the insertion length in CIGAR does not match alt
  res3 <- get_insertion(
    pos      = 15,
    alt      = "ATT",
    r_pos    = 10,
    r_cigar  = "3S6M1I2M",
    r_query  = "AAAAAAAAATTA",
    r_qual   = "#########+##"
  )

  expect_equal(res3$base, "no_insertion_detected")
  expect_equal(res3$qual, "no_insertion_detected")


  # example 4: Case where the position does not correspond to the insertion area
  res4 <- get_insertion(
    pos      = 15,
    alt      = "ATT",
    r_pos    = 9,
    r_cigar  = "3S6M2I1M",
    r_query  = "AAAAAAAAATTA",
    r_qual   = "#########++#"
  )

  expect_equal(res4$base, "no_insertion_detected")
  expect_equal(res4$qual, "no_insertion_detected")

  # example 5: Case where the read cover the position at the last nucleotide (nucleotide before the insertion)
  res5 <- get_insertion(
    pos      = 15,
    alt      = "ATT",
    r_pos    = 9,
    r_cigar  = "3S6M2I1M",
    r_query  = "AAAAAAAAATTA",
    r_qual   = "#########++#"
  )

  expect_equal(res5$base, "no_insertion_detected")
  expect_equal(res5$qual, "no_insertion_detected")

  # example 6: Case with several insertions and the 2nd one good
  res1 <- get_insertion(
    pos      = 15,
    alt      = "ATT",
    r_pos    = 10,
    r_cigar  = "2M2I4M2I1M",
    r_query  = "AAAAAAAATTA",
    r_qual   = "#######I++#"
  )

  expect_equal(res1$base, "+TT")
  expect_equal(res1$qual, "I")
})
