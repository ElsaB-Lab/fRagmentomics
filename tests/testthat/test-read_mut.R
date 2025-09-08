test_that("read mut", {
  vcf_test <- system.file("testdata/mutations/", "mutations_test.vcf", package = "fRagmentomics")
  tsv_test <- system.file("testdata/mutations/", "mutations_test.tsv", package = "fRagmentomics")
  vcf_test_compressed <- system.file("testdata/mutations/", "mutations_test.vcf.gz", package = "fRagmentomics")
  tsv_test_compressed <- system.file("testdata/mutations/", "mutations_test.tsv.gz", package = "fRagmentomics")

  #---------------------------------------
  # chr:pos:ref:alt cases
  #---------------------------------------
  # Normal cases
  expect_equal(read_mut("chr1:123:A:T"), data.frame(CHROM = "chr1", POS = 123, REF = "A", ALT = "T", stringsAsFactors = FALSE))
  # X for Chromosome
  expect_equal(read_mut("chrX:456:G:C"), data.frame(CHROM = "chrX", POS = 456, REF = "G", ALT = "C", stringsAsFactors = FALSE))

  # Special caracters for REF and ALT
  expect_equal(read_mut("4:111:_:-"), data.frame(CHROM = "4", POS = 111, REF = "_", ALT = "-", stringsAsFactors = FALSE))

  # NA in REF of ALT
  expect_equal(read_mut("2:678:NA:T"), data.frame(CHROM = "2", POS = 678, REF = "NA", ALT = "T", stringsAsFactors = FALSE))

  # Mutliallelic cases
  expect_equal(
    read_mut("8:555:A:T,-,"),
    data.frame(
      CHROM = rep("8", 3),
      POS = rep(555, 3),
      REF = rep("A", 3),
      ALT = c("T", "-", "-"),
      stringsAsFactors = FALSE
    )
  )

  # Multiallelic ref is not accepted
  expect_warning(
    expect_error(
      read_mut("8:555:TAC,_,,T:TA"),
      "After reading and filtering, the mutation information is empty. No valid mutation data found."
    ),
    "REF can not be multiallelic"
  )

  # Invalid case to see eror
  # Miss Ref or Alt
  expect_error(
    read_mut("chr1:123:A"),
    "The parameter 'mut' ('chr1:123:A') is not in the expected format (.tsv, .vcf, chr:pos:ref:alt).",
    fixed = TRUE
  )
  # Add extra parameter
  expect_error(
    read_mut("chr1:123:A:T:extra"),
    "The parameter 'mut' ('chr1:123:A:T:extra') is not in the expected format (.tsv, .vcf, chr:pos:ref:alt).",
    fixed = TRUE
  )

  # Invalid_format
  expect_error(
    read_mut("invalid_format"),
    "The parameter 'mut' ('invalid_format') is not in the expected format (.tsv, .vcf, chr:pos:ref:alt).",
    fixed = TRUE
  )

  #---------------------------------------
  # vcf cases
  #---------------------------------------
  # Case 1: VCF
  # Capture the warning to not be shown in the console
  captured_warnings <- testthat::capture_warnings(
    df_vcf_test <- read_mut(vcf_test)
  )

  expected_results <- data.frame(
    CHROM = c(
      "chr1", "2", "chr3", "chr22", "22", "chrX", "chrY", "chr4", "chr5", "chr6",
      "chr8", "chr8", "chrY", "42453", "."
    ),
    POS = c(
      12345, 67890, 101112, 56789, 56789, 54321, 99999, 1234, 5678, 91011,
      151617, 151617, NA, 181920, 181920
    ),
    REF = c(
      "A", "GT", "A.", "NA", ".", "A", ".", "-", ".", "G",
      "G", "G", "T", "AC", "AC"
    ),
    ALT = c(
      "T", "C", "A-", "NA", "-", ".", "C", "A", ".", "NA",
      "A", "T", "C", "-", "-"
    ),
    stringsAsFactors = FALSE
  )

  # Check if the warning ae here
  print(captured_warnings)
  expect_true(any(grepl("non-integer values in the POS column", captured_warnings)))
  expect_true(any(grepl("REF can not be multiallelic", captured_warnings)))
  # Check if the test pass
  expect_equal(df_vcf_test, expected_results)

  # Case 2: Compressed VCF
  # Capture the warning to not be shown in the console
  captured_warnings <- testthat::capture_warnings(
    df_vcf_test2 <- read_mut(vcf_test_compressed)
  )

  expected_results2 <- data.frame(
    CHROM = c(
      "chr1", "2", "chr3", "chr22", "22", "chrX", "chrY", "chr4", "chr5", "chr6",
      "chr8", "chr8", "chrY", "42453", "."
    ),
    POS = c(
      12345, 67890, 101112, 56789, 56789, 54321, 99999, 1234, 5678, 91011,
      151617, 151617, NA, 181920, 181920
    ),
    REF = c(
      "A", "GT", "A.", "NA", ".", "A", ".", "-", ".", "G",
      "G", "G", "T", "AC", "AC"
    ),
    ALT = c(
      "T", "C", "A-", "NA", "-", ".", "C", "A", ".", "NA",
      "A", "T", "C", "-", "-"
    ),
    stringsAsFactors = FALSE
  )

  # Check if the warning ae here
  expect_true(any(grepl("non-integer values in the POS column", captured_warnings)))
  expect_true(any(grepl("REF can not be multiallelic", captured_warnings)))
  # Check if the test pass
  expect_equal(df_vcf_test2, expected_results2)

  #---------------------------------------
  # tsv cases
  #---------------------------------------
  # Case 1: TSV
  captured_warnings <- testthat::capture_warnings(
    df_tsv_test <- read_mut(tsv_test)
  )

  expected_results3 <- data.frame(
    CHROM = c(
      "chr1", "2", "chr3", "chr22", "22", "chrX", "chrY", "chr4", "chr5", "chr6",
      "chr8", "chr8", "chrX", "chrY", "42453", ""
    ),
    POS = c(
      12345, 67890, 101112, 56789, 56789, 54321, 99999, 1234, 5678, 91011,
      151617, 151617, NA, NA, 181920, 181920
    ),
    REF = c(
      "A", "GT", "A.", "NA", ".", "A", "", "-", "", "G",
      "G", "G", "A", "T", "AC", "AC"
    ),
    ALT = c(
      "T", "C", "A-", "NA", "-", "-", "C", "A", "-", "NA",
      "A", "T", "-", "C", "-", "-"
    ),
    stringsAsFactors = FALSE
  )

  # Check if the warning ae here
  expect_true(any(grepl("non-integer values in the POS column", captured_warnings)))
  expect_true(any(grepl("REF can not be multiallelic", captured_warnings)))
  # Check if the test pass
  expect_equal(df_tsv_test, expected_results3)

  # Case 2: Compressed TSV
  captured_warnings <- testthat::capture_warnings(
    df_tsv_test2 <- read_mut(tsv_test_compressed)
  )

  expected_results4 <- data.frame(
    CHROM = c(
      "chr1", "2", "chr3", "chr22", "22", "chrX", "chrY", "chr4", "chr5", "chr6",
      "chr8", "chr8", "chrX", "chrY", "42453", ""
    ),
    POS = c(
      12345, 67890, 101112, 56789, 56789, 54321, 99999, 1234, 5678, 91011,
      151617, 151617, NA, NA, 181920, 181920
    ),
    REF = c(
      "A", "GT", "A.", "NA", ".", "A", "", "-", "", "G",
      "G", "G", "A", "T", "AC", "AC"
    ),
    ALT = c(
      "T", "C", "A-", "NA", "-", "-", "C", "A", "-", "NA",
      "A", "T", "-", "C", "-", "-"
    ),
    stringsAsFactors = FALSE
  )

  # Check if the warning ae here
  expect_true(any(grepl("non-integer values in the POS column", captured_warnings)))
  expect_true(any(grepl("REF can not be multiallelic", captured_warnings)))
  # Check if the test pass
  expect_equal(df_tsv_test2, expected_results4)
})
