test_that("process_fragment_status - SNV with both reads covering", {
  ref <- "C"
  alt <- "A"

  # Example 1: Both reads have the alternative allele -> MUT
  res1 <- process_fragment_status(ref, alt, "SNV", "A", "A")
  expect_equal(res1, "MUT")

  # Example 2: Both reads have the reference allele -> WT
  res2 <- process_fragment_status(ref, alt, "SNV", "C", "C")
  expect_equal(res2, "WT")

  # Example 3: Only one read has the alternative allele -> WT_but_other_read_mut
  res3 <- process_fragment_status(ref, alt, "SNV", "A", "C")
  expect_equal(res3, "WT_but_other_read_mut")

  # Example 4: One read is ref the other is not ref and alt
  res4 <- process_fragment_status(ref, alt, "SNV", "T", "C")
  expect_equal(res4, "WT_but_other_read_mut_with_other_alt")

  # Example 5: One MUT and the other not ref and alt
  res5 <- process_fragment_status(ref, alt, "SNV", "A", "T")
  expect_equal(res5, "MUT_but_other_read_mut_with_other_alt")

  # Example 6: Both read mut with another alt
  res6 <- process_fragment_status(ref, alt, "SNV", "T", "T")
  expect_equal(res6, "Error_both_read_mut_with_other_alt")
})

test_that("process_fragment_status - SNV with incomplete coverage", {
  ref <- "C"
  alt <- "A"

  # Example 1: One read available with the alternative allele -> MUT
  res1 <- process_fragment_status(ref, alt, "SNV", "A", NA)
  expect_equal(res1, "MUT")

  # Example 2: One read available without the alternative allele -> WT
  res2 <- process_fragment_status(ref, alt, "SNV", NA, "C")
  expect_equal(res2, "WT")

  # Example 3: One read available diff than ref and alt
  res3 <- process_fragment_status(ref, alt, "SNV", NA, "T")
  expect_equal(res3, "Other_MUT")

  # Example 4: Both reads missing -> ERROR
  res4 <- process_fragment_status(ref, alt, "SNV", NA, NA)
  expect_equal(res4, "Error_1_read_should_cover_the_position")
})

test_that("process_fragment_status - Insertion with both reads covering", {
  # For insertion type, target value is "insertion_detected"
  ref <- "A"
  alt <- "ATT"

  # Example 1: Both reads detect insertion -> MUT
  res1 <- process_fragment_status(ref, alt, "insertion", "+TT", "+TT")
  expect_equal(res1, "MUT")

  # Example 2: Only one read detects insertion -> WT_but_other_read_mut
  res2 <- process_fragment_status(ref, alt, "insertion", "+TT", "no_insertion_detected")
  expect_equal(res2, "WT_but_other_read_mut")

  # Example 3: Neither read detects insertion -> WT
  res3 <- process_fragment_status(ref, alt, "insertion", "no_insertion_detected", "no_insertion_detected")
  expect_equal(res3, "WT")

  # Example 4: Not supposed to happen because detector-insertion.R only detect the good insertion
  res4 <- process_fragment_status(ref, alt, "insertion", "+TT", "+CC")
  expect_equal(res4, "Error_to_detect_ref/alt")
})

test_that("process_fragment_status - Insertion with incomplete coverage", {
  ref <- "A"
  alt <- "ATT"
  # Example 1: One read available and detects insertion -> MUT
  res1 <- process_fragment_status(ref, alt, "insertion", "+TT", NA)
  expect_equal(res1, "MUT")

  # Example 2: One read available and does not detect insertion -> WT
  res2 <- process_fragment_status(ref, alt, "insertion", NA, "no_insertion_detected")
  expect_equal(res2, "WT")

  # Example 3: Both reads missing -> WT
  res3 <- process_fragment_status(ref, alt, "insertion", NA, NA)
  expect_equal(res3, "Error_1_read_should_cover_the_position")
})

test_that("process_fragment_status - Deletion with both reads covering", {
  ref <- "ATT"
  alt <- "A"

  # Example 1: Both reads detect deletion -> MUT
  res1 <- process_fragment_status(ref, alt, "deletion", "-TT", "-TT")
  expect_equal(res1, "MUT")

  # Example 2: Only one read detects deletion -> WT_but_other_read_mut
  res2 <- process_fragment_status(ref, alt, "deletion", "-TT", "no_deletion_detected")
  expect_equal(res2, "WT_but_other_read_mut")

  # Example 3: Neither read detects deletion -> WT
  res3 <- process_fragment_status(ref, alt, "deletion", "no_deletion_detected", "no_deletion_detected")
  expect_equal(res3, "WT")
})

test_that("process_fragment_status - Deletion with incomplete coverage", {
  ref <- "ATT"
  alt <- "A"

  # Example 1: One read available and detects deletion -> MUT
  res1 <- process_fragment_status(ref, alt, "deletion", "-TT", NA)
  expect_equal(res1, "MUT")

  # Example 2: One read available and does not detect deletion -> WT
  res2 <- process_fragment_status(ref, alt, "deletion", NA, "no_deletion_detected")
  expect_equal(res2, "WT")

  # Example 3: Both reads missing -> WT
  res3 <- process_fragment_status(ref, alt, "deletion", NA, NA)
  expect_equal(res3, "Error_1_read_should_cover_the_position")
})
