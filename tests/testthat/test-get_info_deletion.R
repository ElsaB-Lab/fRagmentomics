test_that("get_info_deletion", {
  # Example 1: Deletion is correctly detected (and not ambiguous by length)
  res1 <- get_info_deletion(
    pos = 9,
    ref = "TAT",
    r_pos = 1,
    r_cigar = "9M2D10M", # Covers 1+9+2+10-1=21. Del at 9+1=10.
    r_qual = "########I##########J", # Qual "I" is at 9th pos
    pos_after_indel_repetition = 12 # threshold = 12-2=10. 21 >= 10. Not ambiguous.
  )
  expect_equal(res1$base, "-AT")
  expect_equal(res1$qual, "I")

  # Example 2: Case where the read doesn't carry the target deletion (not ambiguous by length)
  res2 <- get_info_deletion(
    pos = 9,
    ref = "TA", # Deletion of "A", length 1
    r_pos = 1,
    r_cigar = "15M1D4M", # D op at 15+1=16, not at 9+1=10. Covers 1+15+1+4-1=20.
    r_qual = "####################",
    pos_after_indel_repetition = 11 # threshold = 11-1=10. 20 >= 10. Not ambiguous.
  )
  expect_equal(res2$base, "no_deletion_detected")
  expect_equal(res2$qual, "no_deletion_detected")

  # Example 3: Case where the deletion length in CIGAR does not match (not ambiguous by length)
  res3 <- get_info_deletion(
    pos = 9,
    ref = "TAT", # Expected del "AT", length 2
    r_pos = 1,
    r_cigar = "9M1D5M", # CIGAR has 1D, not 2D. Covers 1+9+1+5-1=15.
    r_qual = "########I#####",
    pos_after_indel_repetition = 12 # threshold = 12-2=10. 15 >= 10. Not ambiguous.
  )
  expect_equal(res3$base, "no_deletion_detected")
  expect_equal(res3$qual, "no_deletion_detected")

  # Example 4: Case where the read doesn't cover the position 'pos'
  res4 <- get_info_deletion(
    pos = 9,
    ref = "TAT",
    r_pos = 1,
    r_cigar = "8M", # Covers 1+8-1=8. pos is 9. Fails coverage check.
    r_qual = "########",
    pos_after_indel_repetition = 11 # Value doesn't matter here
  )
  expect_equal(res4$base, NA_character_)
  expect_equal(res4$qual, NA_character_)

  # Example 5: Read covers up to 'pos', no deletion op, but now can be ambiguous or not.
  # Test 5a: Becomes ambiguous if threshold is higher
  # Ref seq could be: XXXXXXXXTATAC
  res5a <- get_info_deletion(
    pos = 9,
    ref = "TAT", # Deletion of "AT", length 2
    r_pos = 1,
    r_cigar = "9M", # Covers 1+9-1=9. last_ref_pos_covered_by_read = 9.
    r_qual = "########I",
    pos_after_indel_repetition = 12 # threshold = 12-2=10. 9 < 10. Ambiguous.
  )
  expect_equal(res5a$base, "ambiguous")
  expect_equal(res5a$qual, "ambiguous")

  # Test 5b: Stays "no_deletion_detected" if threshold allows
  res5b <- get_info_deletion(
    pos = 9,
    ref = "TAT", # Deletion of "AT", length 2
    r_pos = 1,
    r_cigar = "9M", # Covers 1+9-1=9.
    r_qual = "########I",
    pos_after_indel_repetition = 11 # threshold = 11-2=9. 9 >= 9. Not ambiguous.
  )
  expect_equal(res5b$base, "no_deletion_detected")
  expect_equal(res5b$qual, "no_deletion_detected")

  # Example 6: Case where there are 2 deletions and the 2nd one is correct (not ambiguous by length)
  res6 <- get_info_deletion(
    pos = 9,
    ref = "TAT", # Deletion of "AT", length 2
    r_pos = 1,
    r_cigar = "1M2D6M2D10M", # Del at pos 9. Covers 1+1+2+6+2+10-1=21.
    r_qual = "ABCDEFGHIPQRSTU",
    pos_after_indel_repetition = 12 # threshold = 12-2=10. 21 >= 10. Not ambiguous.
  )
  expect_equal(res6$base, "-AT")
  expect_equal(res6$qual, "G")


  # Ambiguity logic tests
  generic_qual <- paste0(rep("F", 50), collapse = "")

  # N1: Ambiguous - Deletion op PRESENT, but read too short and finish by D
  # Ref: CGACCT / Alt: CG   T
  res_n1 <- get_info_deletion(
    pos = 2, ref = "GACC",
    r_pos = 1,
    r_cigar = "2M3D", # del "ACC", len 3. last_cov=1+2+3+1-1=6
    r_qual = generic_qual,
    pos_after_indel_repetition = 6 # Position of the T - 3
  )
  # Last read covered = 5 - 3 (because D is the last operation) = 2
  # Treshold_ref_pos = 6-3 = 3
  # Treshold_ref_pos > Last read covered so ambiguous
  expect_equal(res_n1$base, "ambiguous")
  expect_equal(res_n1$qual, "ambiguous")

  # N2: Ambiguous - No deletion op, read too short
  res_n2 <- get_info_deletion(
    pos = 10,
    ref = "GACC",
    r_pos = 1,
    r_cigar = "15M", # No D op. last_cov=1+15-1=15
    r_qual = generic_qual, pos_after_indel_repetition = 20 # threshold = 20-3=17. 15 < 17.
  )
  expect_equal(res_n2$base, "ambiguous")
  expect_equal(res_n2$qual, "ambiguous")

  # N3: Deletion Detected (Not Ambiguous) - Deletion op present, read long enough
  res_n3 <- get_info_deletion(
    pos = 10,
    ref = "GACC",
    r_pos = 1,
    r_cigar = "10M3D5M", # Del "ACC", len 3. last_cov=1+10+3+5-1=18
    r_qual = generic_qual, pos_after_indel_repetition = 20 # threshold = 20-3=17. 18 >= 17.
  )
  expect_equal(res_n3$base, "-ACC")
  expect_equal(res_n3$qual, substr(generic_qual, 10, 10))

  # N4: No Deletion Detected (Not Ambiguous) - No deletion op, read long enough
  res_n4 <- get_info_deletion(
    pos = 10,
    ref = "GACC",
    r_pos = 1,
    r_cigar = "17M", # No D op. last_cov=1+18-1=18
    r_qual = generic_qual, pos_after_indel_repetition = 20 # threshold = 20-3=17. 17 >= 17.
  )
  expect_equal(res_n4$base, "no_deletion_detected")
  expect_equal(res_n4$qual, "no_deletion_detected")


  # Tests for warnings and rule skipping

  # N5: pos_after_indel_repetition is NULL - Ambiguity rule not applied
  expect_warning(
    res_n5 <- get_info_deletion(
      pos = 10,
      ref = "GACC",
      r_pos = 1,
      r_cigar = "10M3D2M", # Del op found.
      r_qual = generic_qual, pos_after_indel_repetition = NULL
    ),
    regexp = "'get_deletion': 'pos_after_indel_repetition' is missing, NA, or not numeric. Ambiguity rule based on read length cannot be applied."
  )
  expect_equal(res_n5$base, "-ACC") # Falls back to CIGAR op finding

  # N6: pos_after_indel_repetition is <= 0 - Ambiguity rule not applied
  expect_warning(
    res_n6 <- get_info_deletion( # Assuming get_info_deletion is the correct function name you are testing
      pos = 10,
      ref = "GACC", # del "ACC", len 3
      r_pos = 1,
      r_cigar = "10M3D2M", # Del op found
      r_qual = generic_qual,
      pos_after_indel_repetition = 0
    ),
    # Corrected regexp with escaped parentheses and periods:
    regexp = "'get_deletion': 'pos_after_indel_repetition' \\(0\\) must be positive for the ambiguity rule\\. Rule cannot be applied meaningfully\\."
  )
  expect_equal(res_n6$base, "-ACC") # Falls back to CIGAR op finding as rule isn't applied

  # N7: deletion_length is <= 0 - Ambiguity rule not applied meaningfully
  expect_warning(
    res_n7 <- get_info_deletion( # Assuming get_info_deletion is the correct function name
      pos = 10,
      ref = "G", # deletion_length = 0
      r_pos = 1,
      r_cigar = "12M", # No D op of len 0 will be found by the CIGAR check
      r_qual = generic_qual,
      pos_after_indel_repetition = 20 # Valid pos_after_indel
    ),
    # Corrected regexp with escaped parentheses and periods:
    regexp = "'get_deletion': 'deletion_length' \\(0\\) must be positive for the ambiguity rule related to deletions\\. Rule cannot be applied in this context\\."
  )
  expect_equal(res_n7$base, "no_deletion_detected") # Rule not applied, no D op found.

  # N8: pos_after_indel_repetition is NA
  expect_warning(
    res_n8 <- get_info_deletion(
      pos = 10,
      ref = "GACC",
      r_pos = 1,
      r_cigar = "15M", # No del op
      r_qual = generic_qual, pos_after_indel_repetition = NA_integer_
    ),
    regexp = "'get_deletion': 'pos_after_indel_repetition' is missing, NA, or not numeric. Ambiguity rule based on read length cannot be applied."
  )
  expect_equal(res_n8$base, "no_deletion_detected") # No D op, rule not applied.

  # N9: Deletion at the very start of the read (qual should be NA)
  expect_warning(
    res_n9 <- get_info_deletion( # Assuming get_info_deletion is the correct function name
      pos = 0,
      ref = "GACC",
      r_pos = 1,
      r_cigar = "3D5M",
      r_qual = "FFFFF",
      pos_after_indel_repetition = 10
    ),
    # Corrected regexp with anchors and escaped special characters:
    regexp = "^'get_deletion': 'pos' \\(0\\) is not a valid 1-based coordinate for the base preceding a deletion\\. Returning NA\\.$"
  )
  expect_equal(res_n9$base, NA_character_)
  expect_equal(res_n9$qual, NA_character_)
})
