test_that("get_mutation_status_of_fragment returns correct detail and simple statuses", {
  # --- ERROR Case (Both NA) ---
  result <- get_mutation_status_of_fragment(NA, NA)
  expect_equal(result$Detail, "ERR")
  expect_equal(result$Simple, "ERR")

  # --- MUT Cases ---
  # Both MUT, simple strings
  result <- get_mutation_status_of_fragment("MUT", "MUT")
  expect_equal(result$Detail, "MUT")
  expect_equal(result$Simple, "MUT")

  # MUT with detail, MUT without detail
  result <- get_mutation_status_of_fragment("MUT:C>T", "MUT")
  expect_equal(result$Detail, "MUT & MUT:C>T")
  expect_equal(result$Simple, "MUT")

  # MUT with detail, MUT with different detail
  result <- get_mutation_status_of_fragment("MUT:C>T", "MUT:G>A")
  expect_equal(result$Detail, "MUT:C>T & MUT:G>A")
  expect_equal(result$Simple, "MUT")

  # MUT and MUT potentially larger MUT
  result <- get_mutation_status_of_fragment("MUT", "MUT but potentially larger MUT")
  expect_equal(result$Detail, "MUT & MUT but potentially larger MUT")
  expect_equal(result$Simple, "MUT")

  # AMB and MUT
  result <- get_mutation_status_of_fragment("AMB", "MUT")
  expect_equal(result$Detail, "AMB & MUT")
  expect_equal(result$Simple, "MUT")

  # MUT and AMB (symmetrical check)
  result <- get_mutation_status_of_fragment("MUT", "AMB")
  expect_equal(result$Detail, "AMB & MUT")
  expect_equal(result$Simple, "MUT")

  # NA and MUT
  result <- get_mutation_status_of_fragment(NA, "MUT")
  expect_equal(result$Detail, "MUT")
  expect_equal(result$Simple, "MUT")

  # MUT and NA (symmetrical check)
  result <- get_mutation_status_of_fragment("MUT", NA)
  expect_equal(result$Detail, "MUT")
  expect_equal(result$Simple, "MUT")

  # --- WT Cases ---
  # Both WT, simple strings
  result <- get_mutation_status_of_fragment("WT", "WT")
  expect_equal(result$Detail, "WT")
  expect_equal(result$Simple, "WT")

  # WT with detail, WT without detail
  result <- get_mutation_status_of_fragment("WT:ref", "WT")
  expect_equal(result$Detail, "WT & WT:ref")
  expect_equal(result$Simple, "WT")

  # AMB and WT
  result <- get_mutation_status_of_fragment("AMB", "WT")
  expect_equal(result$Detail, "AMB & WT")
  expect_equal(result$Simple, "WT")

  # NA and WT
  result <- get_mutation_status_of_fragment(NA, "WT")
  expect_equal(result$Detail, "WT")
  expect_equal(result$Simple, "WT")

  # --- OTH Cases ---
  # Both OTH, simple strings
  result <- get_mutation_status_of_fragment("OTH", "OTH")
  expect_equal(result$Detail, "OTH")
  expect_equal(result$Simple, "OTH")

  # OTH with detail, OTH without detail
  result <- get_mutation_status_of_fragment("OTH:indel", "OTH")
  expect_equal(result$Detail, "OTH & OTH:indel")
  expect_equal(result$Simple, "OTH")

  # AMB and OTH
  result <- get_mutation_status_of_fragment("AMB", "OTH")
  expect_equal(result$Detail, "AMB & OTH")
  expect_equal(result$Simple, "OTH")

  # NA and OTH
  result <- get_mutation_status_of_fragment(NA, "OTH")
  expect_equal(result$Detail, "OTH")
  expect_equal(result$Simple, "OTH")

  # --- Discordant Cases ---
  # MUT and WT
  result <- get_mutation_status_of_fragment("MUT", "WT")
  expect_equal(result$Detail, "MUT & WT")
  expect_equal(result$Simple, "N/I")

  # WT and MUT (symmetrical check)
  result <- get_mutation_status_of_fragment("WT", "MUT")
  expect_equal(result$Detail, "MUT & WT")
  expect_equal(result$Simple, "N/I")

  # MUT with detail and WT with detail
  result <- get_mutation_status_of_fragment("MUT:C>T", "WT:ref")
  expect_equal(result$Detail, "MUT:C>T & WT:ref")
  expect_equal(result$Simple, "N/I")

  # MUT and OTH
  result <- get_mutation_status_of_fragment("MUT", "OTH")
  expect_equal(result$Detail, "MUT & OTH")
  expect_equal(result$Simple, "N/I")

  # OTH and MUT (symmetrical check)
  result <- get_mutation_status_of_fragment("OTH", "MUT")
  expect_equal(result$Detail, "MUT & OTH")
  expect_equal(result$Simple, "N/I")

  # WT and OTH
  result <- get_mutation_status_of_fragment("WT", "OTH")
  expect_equal(result$Detail, "OTH & WT")
  expect_equal(result$Simple, "N/I")

  # OTH and WT (symmetrical check)
  result <- get_mutation_status_of_fragment("OTH", "WT")
  expect_equal(result$Detail, "OTH & WT")
  expect_equal(result$Simple, "N/I")

  # --- Discordant Cases "potentially" ---
  # WT & uncertainty for MUT -> WT
  result <- get_mutation_status_of_fragment("WT", "MUT by CIGAR but potentially WT")
  expect_equal(result$Detail, "MUT by CIGAR but potentially WT & WT")
  expect_equal(result$Simple, "WT")

  #  MUT & uncertainty for WT -> MUT
  result <- get_mutation_status_of_fragment("MUT", "WT by CIGAR but potentially MUT")
  expect_equal(result$Detail, "MUT & WT by CIGAR but potentially MUT")
  expect_equal(result$Simple, "MUT")

  # MUT & uncertainty for OTH and WT -> N/I (no proof towards MUT)
  result <- get_mutation_status_of_fragment("MUT", "WT by CIGAR but potentially OTH")
  expect_equal(result$Detail, "MUT & WT by CIGAR but potentially OTH")
  expect_equal(result$Simple, "N/I")

  # MUT & WT but both uncertainty for OTH -> N/I (both uncertainty, no one to trust)
  result <- get_mutation_status_of_fragment("MUT by CIGAR but potentially OTH", "WT by CIGAR but potentially OTH")
  expect_equal(result$Detail, "MUT by CIGAR but potentially OTH & WT by CIGAR but potentially OTH")
  expect_equal(result$Simple, "N/I")

  # --- Ambiguous Cases ---
  # AMB and AMB
  result <- get_mutation_status_of_fragment("AMB", "AMB")
  expect_equal(result$Detail, "AMB")
  expect_equal(result$Simple, "N/I")

  # AMB with detail and AMB without detail
  result <- get_mutation_status_of_fragment("AMB:low_qual", "AMB")
  expect_equal(result$Detail, "AMB & AMB:low_qual")
  expect_equal(result$Simple, "N/I")
})
