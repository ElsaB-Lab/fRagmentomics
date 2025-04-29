test_that("get deletion", {
  # example 1 : deletion is correctly detected
  res1 <- get_deletion(
    pos     = 9,
    ref     = "TAT",
    r_pos   = 1,
    r_cigar = "9M2D10M",
    r_qual  = "########I##########"
  )

  expect_equal(res1$base, "-AT")
  expect_equal(res1$qual, "I")

  # example 2: Case where the read doesn't carry the deletion
  res2 <- get_deletion(
    pos     = 9,
    ref     = "TA",
    r_pos   = 1,
    r_cigar = "15M1D4M",
    r_qual  = "########I##########"
  )

  expect_equal(res2$base, "no_deletion_detected")
  expect_equal(res2$qual, "no_deletion_detected")

  # example 3: Case where the deletion length does not match the expected length
  res3 <- get_deletion(
    pos     = 9,
    ref     = "TAT",
    r_pos   = 1,
    r_cigar = "8M2D",
    r_qual  = "########I##########"
  )

  expect_equal(res3$base, "no_deletion_detected")
  expect_equal(res3$qual, "no_deletion_detected")

  # example 4: Case where the read doesn't cover the position of the deletion
  res4 <- get_deletion(
    pos     = 9,
    ref     = "TAT",
    r_pos   = 1,
    r_cigar = "8M",
    r_qual  = "########I##########"
  )

  expect_equal(res4$base, NA)
  expect_equal(res4$qual, NA)

  # example 5: Case where the read cover the position at the last nucleotide (nucleotide before the deletion)
  res5 <- get_deletion(
    pos     = 9,
    ref     = "TAT",
    r_pos   = 1,
    r_cigar = "9M",
    r_qual  = "########I##########"
  )

  expect_equal(res5$base, "no_deletion_detected")
  expect_equal(res5$qual, "no_deletion_detected")

  # example 6: Case where there are 2 deletions and the 2nd one is correct
  res1 <- get_deletion(
    pos     = 9,
    ref     = "TAT",
    r_pos   = 1,
    r_cigar = "1M2D6M2D10M",
    r_qual  = "######I##########"
  )

  expect_equal(res1$base, "-AT")
  expect_equal(res1$qual, "I")
})
