# ======================================================================
# Tests for SNVs (Single Nucleotide Variants)
# ======================================================================

test_that("process_fragment_status - SNV both covering", {
  ref <- "C"
  alt <- "A"

  # MUT + MUT -> MUT
  expect_equal(process_fragment_status(ref, alt, "SNV", "A", "A"), "MUT")

  # WT + WT -> WT
  expect_equal(process_fragment_status(ref, alt, "SNV", "C", "C"), "WT")

  # MUT + WT -> WT_but_other_read_mut
  expect_equal(process_fragment_status(ref, alt, "SNV", "A", "C"), "WT_but_other_read_mut")
  expect_equal(process_fragment_status(ref, alt, "SNV", "C", "A"), "WT_but_other_read_mut")

  # WT + other -> WT_but_other_read_mut_with_other_alt
  expect_equal(process_fragment_status(ref, alt, "SNV", "T", "C"), "WT_but_other_read_mut_with_other_alt")

  # MUT + other -> MUT_but_other_read_mut_with_other_alt
  expect_equal(process_fragment_status(ref, alt, "SNV", "A", "T"), "MUT_but_other_read_mut_with_other_alt")

  # other + other -> Error_both_read_mut_with_other_alt
  expect_equal(process_fragment_status(ref, alt, "SNV", "T", "G"), "Error_both_read_mut_with_other_alt")
})

test_that("process_fragment_status - SNV one read covering", {
  ref <- "C"
  alt <- "A"

  # MUT + NA -> MUT
  expect_equal(process_fragment_status(ref, alt, "SNV", "A", NA), "MUT")

  # WT + NA -> WT
  expect_equal(process_fragment_status(ref, alt, "SNV", NA, "C"), "WT")

  # other + NA -> Other_MUT
  expect_equal(process_fragment_status(ref, alt, "SNV", NA, "T"), "Other_MUT")

  # NA + NA -> Error
  expect_equal(process_fragment_status(ref, alt, "SNV", NA, NA), "Error_1_read_should_cover_the_position")
})

# ======================================================================
# Tests for Indels
# ======================================================================

test_that("process_fragment_status - Indels", {
  # --- Insertion ---
  ref_ins <- "A"
  alt_ins <- "ATT" # Mutation = +TT

  # MUT + MUT -> MUT
  expect_equal(process_fragment_status(ref_ins, alt_ins, "insertion", "+TT", "+TT"), "MUT")
  # MUT + WT -> WT_but_other_read_mut
  expect_equal(process_fragment_status(ref_ins, alt_ins, "insertion", "+TT", "no_insertion_detected"), "WT_but_other_read_mut")
  # WT + WT -> WT
  expect_equal(process_fragment_status(ref_ins, alt_ins, "insertion", "no_insertion_detected", "no_insertion_detected"), "WT")
  # MUT + NA -> MUT
  expect_equal(process_fragment_status(ref_ins, alt_ins, "insertion", "+TT", NA), "MUT")
  # WT + NA -> WT
  expect_equal(process_fragment_status(ref_ins, alt_ins, "insertion", NA, "no_insertion_detected"), "WT")

  # --- Deletion ---
  ref_del <- "ATT"
  alt_del <- "A" # Mutation = -TT

  # MUT + MUT -> MUT
  expect_equal(process_fragment_status(ref_del, alt_del, "deletion", "-TT", "-TT"), "MUT")
  # MUT + WT -> WT_but_other_read_mut
  expect_equal(process_fragment_status(ref_del, alt_del, "deletion", "-TT", "no_deletion_detected"), "WT_but_other_read_mut")
  # WT + WT -> WT
  expect_equal(process_fragment_status(ref_del, alt_del, "deletion", "no_deletion_detected", "no_deletion_detected"), "WT")
  # MUT + NA -> MUT
  expect_equal(process_fragment_status(ref_del, alt_del, "deletion", "-TT", NA), "MUT")
  # WT + NA -> WT
  expect_equal(process_fragment_status(ref_del, alt_del, "deletion", NA, "no_deletion_detected"), "WT")
})

test_that("process_fragment_status - Indels with 'ambiguous'", {
  # --- Insertion ---
  ref_ins <- "A"
  alt_ins <- "ATT" # Mutation = +TT

  # Ambiguous + Ambiguous -> Ambiguous
  expect_equal(process_fragment_status(ref_ins, alt_ins, "insertion", "ambiguous", "ambiguous"), "Ambiguous")
  # Ambiguous + NA -> Ambiguous
  expect_equal(process_fragment_status(ref_ins, alt_ins, "insertion", "ambiguous", NA), "Ambiguous")
  # WT + Ambiguous -> WT
  expect_equal(process_fragment_status(ref_ins, alt_ins, "insertion", "no_insertion_detected", "ambiguous"), "WT")
  # MUT + Ambiguous -> MUT
  expect_equal(process_fragment_status(ref_ins, alt_ins, "insertion", "+TT", "ambiguous"), "MUT")

  # --- Deletion ---
  ref_del <- "ATT"
  alt_del <- "A" # Mutation = -TT

  # Ambiguous + Ambiguous -> Ambiguous
  expect_equal(process_fragment_status(ref_del, alt_del, "deletion", "ambiguous", "ambiguous"), "Ambiguous")
  # Ambiguous + NA -> Ambiguous
  expect_equal(process_fragment_status(ref_del, alt_del, "deletion", NA, "ambiguous"), "Ambiguous")
  # WT + Ambiguous -> WT
  expect_equal(process_fragment_status(ref_del, alt_del, "deletion", "ambiguous", "no_deletion_detected"), "WT")
  # MUT + Ambiguous -> MUT
  expect_equal(process_fragment_status(ref_del, alt_del, "deletion", "ambiguous", "-TT"), "MUT")
})


test_that("process_fragment_status - Indels with 'other'", {
  # --- Insertion ---
  ref_ins <- "A"
  alt_ins <- "ATT" # Mutation = +TT

  # WT + other -> WT_but_other_read_mut_with_other_alt
  expect_equal(process_fragment_status(ref_ins, alt_ins, "insertion", "no_insertion_detected", "+GG"), "WT_but_other_read_mut_with_other_alt")

  # MUT + other -> MUT_but_other_read_mut_with_other_alt
  expect_equal(process_fragment_status(ref_ins, alt_ins, "insertion", "+TT", "+GG"), "MUT_but_other_read_mut_with_other_alt")

  # other + other -> Error_both_read_mut_with_other_alt
  expect_equal(process_fragment_status(ref_ins, alt_ins, "insertion", "+CC", "+GG"), "Error_both_read_mut_with_other_alt")

  # other + NA -> Other_MUT
  expect_equal(process_fragment_status(ref_ins, alt_ins, "insertion", "+GG", NA), "Other_MUT")
})
