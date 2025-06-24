test_that("get_insertion handles VCF-consistent alleles and ambiguity correctly", {
  # --- Test Setup ---
  generic_qual <- paste0(rep("I", 50), collapse = "")

  # --- Group 1: Basic Scenarios (Not Ambiguous) ---

  # Example 1: Insertion is correctly detected
  res1 <- get_insertion(
    pos = 15,
    alt = "CAG", # VCF ALT -> REF base is 'C', pure insertion is "AG" (length 2)
    r_pos = 10,
    r_cigar = "6M2I5M", # CIGAR op is 2I
    r_query = "TTCGCCAGGTACCT", # The 6th base is now 'C', followed by inserted "AG"
    r_qual = "FFFFF#GGGGGGGG", # Quality of 6th base ('C') is '#'
    pos_after_indel_repetition = 20 # last_cov=10+6+5-1=20. Not ambiguous.
  )
  expect_equal(res1$base, "+AG")
  expect_equal(res1$qual, "#")

  # Example 2: Preceding base in read does NOT match 'alt'
  res2 <- get_insertion(
    pos = 15,
    alt = "CAG", # Expects preceding base to be 'C'
    r_pos = 10,
    r_cigar = "6M2I5M",
    r_query = "TTCGCTAGGTACC", # Preceding base is 'T', which is incorrect.
    r_qual = "FFFFF#GGGGGGGG",
    pos_after_indel_repetition = 20
  )
  expect_equal("no_insertion_detected", "no_insertion_detected")

  # Example 3: CIGAR insertion length does not match
  res3 <- get_insertion(
    pos = 15,
    alt = "CAG", # Expected pure insertion "AG", length 2
    r_pos = 10,
    r_cigar = "6M3I5M", # CIGAR has 3I, not 2I
    r_query = "TTCGCCAGGGTACC", # Read has 3 bases inserted
    r_qual = generic_qual,
    pos_after_indel_repetition = 20
  )
  expect_equal(res3$base, "no_insertion_detected")

  # Example 4: Insertion is correctly detected but not the good one
  res4 <- get_insertion(
    pos = 15,
    alt = "CAG", # VCF ALT -> REF base is 'C', pure insertion is "AG" (length 2)
    r_pos = 10,
    r_cigar = "6M2I5M", # CIGAR op is 2I
    r_query = "TTCGCCATGTACCT",
    r_qual = "FFFFF#GGGGGGGG",
    pos_after_indel_repetition = 20 # last_cov=10+6+5-1=20. Not ambiguous.
  )
  expect_equal(res4$base, "+AT")
  expect_equal(res4$qual, "#")

  # Example 5: Position of interest is not covered (too high)
  res5 <- get_insertion(
    pos = 30,
    alt = "CAG",
    r_pos = 10,
    r_cigar = "6M2I5M",
    r_query = "TTCGCCAGGTACCT",
    r_qual = "FFFFF#GGGGGGGG",
    pos_after_indel_repetition = 25
  )
  expect_equal(res5$base, NA_character_)
  expect_equal(res5$qual, NA_character_)

  # Example 6: Position of interest is not covered (too short)
  res6 <- get_insertion(
    pos = 5,
    alt = "CAG",
    r_pos = 10,
    r_cigar = "6M2I5M",
    r_query = "TTCGCCAGGTACCT",
    r_qual = "FFFFF#GGGGGGGG",
    pos_after_indel_repetition = 25
  )
  expect_equal(res6$base, NA_character_)
  expect_equal(res6$qual, NA_character_)

  # --- Group 2: Ambiguity Logic Tests ---

  # N1: Ambiguous - Insertion op PRESENT, but read too short
  res_n1 <- get_insertion(
    pos = 15, alt = "CAG", r_pos = 10, r_cigar = "6M2I5M", # last_cov = 20
    r_query = "TTCGCCAGAGA", r_qual = generic_qual, # Correct read sequence
    pos_after_indel_repetition = 21 # threshold = 21. 20 < 21 -> Ambiguous
  )
  expect_equal(res_n1$base, "ambiguous")
  expect_equal(res_n1$qual, "ambiguous")

  # N2: Ambiguous - No insertion op, read too short
  res_n2 <- get_insertion(
    pos = 15, alt = "CAG", r_pos = 10, r_cigar = "10M", # last_cov = 19
    r_query = "TTCGCCAGAG", r_qual = generic_qual,
    pos_after_indel_repetition = 20 # threshold = 20. 19 < 20 -> Ambiguous
  )
  expect_equal(res_n2$base, "ambiguous")
  expect_equal(res_n2$qual, "ambiguous")

  # --- Group 3: Warnings and Edge Case Tests ---

  # N3: 'alt' does not represent an insertion (length <= 1)
  expect_warning(
    res_n3 <- get_insertion(
      pos = 15, alt = "C", r_pos = 10, r_cigar = "10M", r_query = "TTCGCTGTAC",
      r_qual = generic_qual, pos_after_indel_repetition = 20
    ),
    regexp = "'get_insertion': 'alt' sequence \\(insertion\\) has length 0"
  )
  expect_equal(res_n3$base, "no_insertion_detected")

  # N4: Insertion at the start of the read alignment (preceding base check is skipped)
  res_n4 <- get_insertion(
    pos = 9,
    alt = "TCAG", # pure insertion "CAG", length 3. Ref base is 'T'.
    r_pos = 10,
    r_cigar = "3I10M", # Insertion right after pos 9
    r_query = "CAGTTCGCTGTAC", # Read starts with insertion, no preceding base to check
    r_qual = "###IIIIIIIIII",
    pos_after_indel_repetition = 5
  )
  expect_equal(res_n4$base, NA_character_)
  expect_equal(res_n4$qual, NA_character_)

  # N5: Ambiguous - Insertion op PRESENT, but read too short (limit)
  res_n5 <- get_insertion(
    # ref = TTCGCCT
    pos = 15, alt = "CA", r_pos = 10, r_cigar = "6M", # last_cov = 20
    r_query = "TTCGCC", r_qual = generic_qual, # Correct read sequence
    pos_after_indel_repetition = 16 # T after the insertion
  )
  # last_nucleotide_read < pos_after_indel_repetition
  expect_equal(res_n5$base, "ambiguous")
  expect_equal(res_n5$qual, "ambiguous")
})
