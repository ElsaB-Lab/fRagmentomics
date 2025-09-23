test_that("get_fragment_mutation_statuses returns correct detail and simple statuses", {

  # --- ERROR Case (Both NA) ---
  result <- get_fragment_mutation_statuses(NA, NA)
  expect_equal(result$Detail, "ERR")
  expect_equal(result$Simple, "ERR")

  # --- MUT Cases ---
  # Both MUT, simple strings
  result <- get_fragment_mutation_statuses("MUT", "MUT")
  expect_equal(result$Detail, "MUT")
  expect_equal(result$Simple, "MUT")

  # MUT with detail, MUT without detail
  result <- get_fragment_mutation_statuses("MUT:C>T", "MUT")
  expect_equal(result$Detail, "MUT & MUT:C>T")
  expect_equal(result$Simple, "MUT")

  # MUT with detail, MUT with different detail
  result <- get_fragment_mutation_statuses("MUT:C>T", "MUT:G>A")
  expect_equal(result$Detail, "MUT:C>T & MUT:G>A")
  expect_equal(result$Simple, "MUT")

  # MUT and MUT potentially larger MUT
  result <- get_fragment_mutation_statuses("MUT", "MUT but potentially larger MUT")
  expect_equal(result$Detail, "MUT:C>T & MUT:G>A")
  expect_equal(result$Simple, "MUT")

  # AMB and MUT
  result <- get_fragment_mutation_statuses("AMB", "MUT")
  expect_equal(result$Detail, "AMB & MUT")
  expect_equal(result$Simple, "MUT")

  # MUT and AMB (symmetrical check)
  result <- get_fragment_mutation_statuses("MUT", "AMB")
  expect_equal(result$Detail, "AMB & MUT")
  expect_equal(result$Simple, "MUT")

  # NA and MUT
  result <- get_fragment_mutation_statuses(NA, "MUT")
  expect_equal(result$Detail, "MUT")
  expect_equal(result$Simple, "MUT")

  # MUT and NA (symmetrical check)
  result <- get_fragment_mutation_statuses("MUT", NA)
  expect_equal(result$Detail, "MUT")
  expect_equal(result$Simple, "MUT")

  # --- WT Cases ---
  # Both WT, simple strings
  result <- get_fragment_mutation_statuses("WT", "WT")
  expect_equal(result$Detail, "WT")
  expect_equal(result$Simple, "WT")

  # WT with detail, WT without detail
  result <- get_fragment_mutation_statuses("WT:ref", "WT")
  expect_equal(result$Detail, "WT & WT:ref")
  expect_equal(result$Simple, "WT")

  # AMB and WT
  result <- get_fragment_mutation_statuses("AMB", "WT")
  expect_equal(result$Detail, "AMB & WT")
  expect_equal(result$Simple, "WT")

  # NA and WT
  result <- get_fragment_mutation_statuses(NA, "WT")
  expect_equal(result$Detail, "WT")
  expect_equal(result$Simple, "WT")

  # --- OTH Cases ---
  # Both OTH, simple strings
  result <- get_fragment_mutation_statuses("OTH", "OTH")
  expect_equal(result$Detail, "OTH")
  expect_equal(result$Simple, "OTH")

  # OTH with detail, OTH without detail
  result <- get_fragment_mutation_statuses("OTH:indel", "OTH")
  expect_equal(result$Detail, "OTH & OTH:indel")
  expect_equal(result$Simple, "OTH")

  # AMB and OTH
  result <- get_fragment_mutation_statuses("AMB", "OTH")
  expect_equal(result$Detail, "AMB & OTH")
  expect_equal(result$Simple, "OTH")

  # NA and OTH
  result <- get_fragment_mutation_statuses(NA, "OTH")
  expect_equal(result$Detail, "OTH")
  expect_equal(result$Simple, "OTH")

  # --- Discordant Cases ---
  # MUT and WT
  result <- get_fragment_mutation_statuses("MUT", "WT")
  expect_equal(result$Detail, "MUT & WT")
  expect_equal(result$Simple, "DIS")

  # WT and MUT (symmetrical check)
  result <- get_fragment_mutation_statuses("WT", "MUT")
  expect_equal(result$Detail, "MUT & WT")
  expect_equal(result$Simple, "DIS")

  # MUT with detail and WT with detail
  result <- get_fragment_mutation_statuses("MUT:C>T", "WT:ref")
  expect_equal(result$Detail, "MUT:C>T & WT:ref")
  expect_equal(result$Simple, "DIS")

  # MUT and OTH
  result <- get_fragment_mutation_statuses("MUT", "OTH")
  expect_equal(result$Detail, "MUT & OTH")
  expect_equal(result$Simple, "DIS")

  # OTH and MUT (symmetrical check)
  result <- get_fragment_mutation_statuses("OTH", "MUT")
  expect_equal(result$Detail, "MUT & OTH")
  expect_equal(result$Simple, "DIS")

  # WT and OTH
  result <- get_fragment_mutation_statuses("WT", "OTH")
  expect_equal(result$Detail, "OTH & WT")
  expect_equal(result$Simple, "DIS")

  # OTH and WT (symmetrical check)
  result <- get_fragment_mutation_statuses("OTH", "WT")
  expect_equal(result$Detail, "OTH & WT")
  expect_equal(result$Simple, "DIS")

  # --- Ambiguous Cases ---
  # AMB and AMB
  result <- get_fragment_mutation_statuses("AMB", "AMB")
  expect_equal(result$Detail, "AMB")
  expect_equal(result$Simple, "AMB")

  # AMB with detail and AMB without detail
  result <- get_fragment_mutation_statuses("AMB:low_qual", "AMB")
  expect_equal(result$Detail, "AMB & AMB:low_qual")
  expect_equal(result$Simple, "AMB")

})
