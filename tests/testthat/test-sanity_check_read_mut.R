test_that("sanity check read mut", {
  # --- 1. Test ---
  valid_mut_info <- data.frame(
    CHROM = c(
      "chr1", "2", "chr3", "chr22", "22", "chrX", "chrY", "chr4", "chr5", "chr6",
      "chr7", "chr8", "42453", ""
    ),
    POS = c(
      12345, 67890, 101112, 56789, 56789, 54321, 99999, 1234, 5678, 91011,
      121314, 151617, 181920, 181920
    ),
    REF = c(
      "A", "GT", "A.", "NA", ".", "A", "", "-", "", "G",
      "T,T", "G", "AC", "AC"
    ),
    ALT = c(
      "T", "C", "A-", "NA", "-", "", "C", "A", "", NA,
      "C", "A,T", "-", "-"
    ),
    stringsAsFactors = FALSE
  )

  expected_results <- data.frame(
    CHROM = c("chr1", "2", "chr3", "chrX", "chrY", "chr4", "chr6"),
    POS = c(12345, 67890, 101112, 54321, 99999, 1234, 91011),
    REF = c("A", "GT", "A.", "A", "", "-", "G"),
    ALT = c("T", "C", "A-", "", "C", "A", NA),
    stringsAsFactors = FALSE
  )

  # Warning capture
  warnings <- testthat::capture_warnings({
    actual_results <- sanity_check_read_mut(valid_mut_info)

    # Reinitiate rownames to avoid error
    rownames(actual_results) <- NULL
    rownames(expected_results) <- NULL

    # Check if equal
    expect_equal(actual_results, expected_results)
  })

  # Check if warnings has been generated
  expect_true(any(grepl("Invalid row", warnings)),
    info = "At least one warning 'Invalid row' was expected"
  )

  # --- 2. Test si POS est string ---
  valid_mut_info2 <- data.frame(
    CHROM = c("chr9", "chrX", "chrY"),
    POS = c("181920", "", "."),
    REF = c("AC", "A", "T"),
    ALT = c("-", "-", "C"),
    stringsAsFactors = FALSE
  )

  # Here we are espected a specific error
  expect_error(
    suppressWarnings(sanity_check_read_mut(valid_mut_info2)),
    "No valid mutations found after sanity check."
  )
})
