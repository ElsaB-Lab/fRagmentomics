library(testthat)

test_that("process_fragment_status - SNV with both reads covering", {
  alt <- "A"

  # Example 1: Both reads have the alternative allele -> MUT
  res1 <- process_fragment_status(alt, "SNV", "A", "A")
  expect_equal(res1, "MUT")

  # Example 2: Only one read has the alternative allele -> WT_but_other_read_mut
  res2 <- process_fragment_status(alt, "SNV", "A", "C")
  expect_equal(res2, "WT_but_other_read_mut")

  # Example 3: Neither read has the alternative allele -> WT
  res3 <- process_fragment_status(alt, "SNV", "T", "C")
  expect_equal(res3, "WT")
})

test_that("process_fragment_status - SNV with incomplete coverage", {
  alt <- "A"

  # Example 1: One read available with the alternative allele -> MUT
  res1 <- process_fragment_status(alt, "SNV", "A", NA)
  expect_equal(res1, "MUT")

  # Example 2: One read available without the alternative allele -> WT
  res2 <- process_fragment_status(alt, "SNV", NA, "C")
  expect_equal(res2, "WT")

  # Example 3: Both reads missing -> WT
  res3 <- process_fragment_status(alt, "SNV", NA, NA)
  expect_equal(res3, "WT")
})

test_that("process_fragment_status - Insertion with both reads covering", {
  # For insertion type, target value is "insertion_detected"

  # Example 1: Both reads detect insertion -> MUT
  res1 <- process_fragment_status("A", "insertion", "+TT", "+GG")
  expect_equal(res1, "MUT")

  # Example 2: Only one read detects insertion -> WT_but_other_read_mut
  res2 <- process_fragment_status("A", "insertion", "+TT", "C")
  expect_equal(res2, "WT_but_other_read_mut")

  # Example 3: Neither read detects insertion -> WT
  res3 <- process_fragment_status("A", "insertion", "C", "G")
  expect_equal(res3, "WT")
})

test_that("process_fragment_status - Insertion with incomplete coverage", {
  # Example 1: One read available and detects insertion -> MUT
  res1 <- process_fragment_status("A", "insertion", "+TT", NA)
  expect_equal(res1, "MUT")

  # Example 2: One read available and does not detect insertion -> WT
  res2 <- process_fragment_status("A", "insertion", NA, "C")
  expect_equal(res2, "WT")

  # Example 3: Both reads missing -> WT
  res3 <- process_fragment_status("A", "insertion", NA, NA)
  expect_equal(res3, "WT")
})

test_that("process_fragment_status - Deletion with both reads covering", {
  # Example 1: Both reads detect deletion -> MUT
  res1 <- process_fragment_status("A", "deletion", "-AC", "-T")
  expect_equal(res1, "MUT")

  # Example 2: Only one read detects deletion -> WT_but_other_read_mut
  res2 <- process_fragment_status("A", "deletion", "-G", "G")
  expect_equal(res2, "WT_but_other_read_mut")

  # Example 3: Neither read detects deletion -> WT
  res3 <- process_fragment_status("A", "deletion", "G", "A")
  expect_equal(res3, "WT")
})


test_that("process_fragment_status - Deletion with incomplete coverage", {
  # Example 1: One read available and detects deletion -> MUT
  res1 <- process_fragment_status("A", "deletion", "-AG", NA)
  expect_equal(res1, "MUT")

  # Example 2: One read available and does not detect deletion -> WT
  res2 <- process_fragment_status("A", "deletion", NA, "C")
  expect_equal(res2, "WT")

  # Example 3: Both reads missing -> WT
  res3 <- process_fragment_status("A", "deletion", NA, NA)
  expect_equal(res3, "WT")
})
