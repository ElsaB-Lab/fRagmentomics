test_that("get_fragment_mutation_status_broad works", {
  
  # MUT Case
  fragment_status_broad <- get_fragment_mutation_status_broad ("MUT", "MUT")
  expect_equal(fragment_status_broad, "MUT")

  fragment_status_broad <- get_fragment_mutation_status_broad ("AMB", "MUT")
  expect_equal(fragment_status_broad, "MUT")

  fragment_status_broad <- get_fragment_mutation_status_broad (NA, "MUT")
  expect_equal(fragment_status_broad, "MUT")

  # WT Case
  fragment_status_broad <- get_fragment_mutation_status_broad ("WT", "WT")
  expect_equal(fragment_status_broad, "WT")

  fragment_status_broad <- get_fragment_mutation_status_broad ("AMB", "WT")
  expect_equal(fragment_status_broad, "WT")

  fragment_status_broad <- get_fragment_mutation_status_broad (NA, "WT")
  expect_equal(fragment_status_broad, "WT")

  # Other MUT Case
  fragment_status_broad <- get_fragment_mutation_status_broad ("Other MUT", "Other MUT")
  expect_equal(fragment_status_broad, "Other MUT")

  fragment_status_broad <- get_fragment_mutation_status_broad ("AMB", "Other MUT")
  expect_equal(fragment_status_broad, "Other MUT")

  fragment_status_broad <- get_fragment_mutation_status_broad (NA, "Other MUT")
  expect_equal(fragment_status_broad, "Other MUT")

  # Discordant cases
  fragment_status_broad <- get_fragment_mutation_status_broad ("MUT", "WT")
  expect_equal(fragment_status_broad, "WT & MUT")

  fragment_status_broad <- get_fragment_mutation_status_broad ("WT", "MUT")
  expect_equal(fragment_status_broad, "WT & MUT")

  fragment_status_broad <- get_fragment_mutation_status_broad ("MUT", "Other MUT")
  expect_equal(fragment_status_broad, "Other MUT & MUT")
  
  fragment_status_broad <- get_fragment_mutation_status_broad ("WT", "Other MUT")
  expect_equal(fragment_status_broad, "WT & Other MUT")

  # Ambiguous Case
  fragment_status_broad <- get_fragment_mutation_status_broad ("AMB", NA)
  expect_equal(fragment_status_broad, "AMB")
  
  fragment_status_broad <- get_fragment_mutation_status_broad ("AMB", "AMB")
  expect_equal(fragment_status_broad, "AMB")

  # NA Case
  fragment_status_broad <- get_fragment_mutation_status_broad (NA, NA)
  expect_equal(fragment_status_broad, "ERROR: position not covered")
})
