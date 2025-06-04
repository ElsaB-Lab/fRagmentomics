test_that("get_fragment_mutation_status works", {
  
  # Non-target MUT Case
  fragment_status <- get_fragment_mutation_status ("WT")
  expect_equal(fragment_status, "Non-target MUT")

  fragment_status <- get_fragment_mutation_status ("Other MUT")
  expect_equal(fragment_status, "Non-target MUT")

  fragment_status <- get_fragment_mutation_status ("WT & Other MUT")
  expect_equal(fragment_status, "Non-target MUT")

  # MUT Case
  fragment_status <- get_fragment_mutation_status ("MUT")
  expect_equal(fragment_status, "MUT")

  # DISCORDANT Case
  fragment_status <- get_fragment_mutation_status ("WT & MUT")
  expect_equal(fragment_status, "DISCORDANT")

  fragment_status <- get_fragment_mutation_status ("Other MUT & MUT")
  expect_equal(fragment_status, "DISCORDANT")

  # AMB cases
  fragment_status <- get_fragment_mutation_status ("AMB")
  expect_equal(fragment_status, "AMB")

  # Error cases
  fragment_status <- get_fragment_mutation_status ("ERROR: position not covered")
  expect_equal(fragment_status, "ERROR: position not covered")
})