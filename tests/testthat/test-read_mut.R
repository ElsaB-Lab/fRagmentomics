test_that("read_mut works with character string input", {
  #---------------------------------------
  # chr:pos:ref:alt cases
  #---------------------------------------
  # Normal cases
  expect_equal(read_mut("chr1:123:A:T"), data.frame(CHROM = "chr1", POS = 123, REF = "A", ALT = "T", stringsAsFactors = FALSE))
  # X for Chromosome
  expect_equal(read_mut("chrX:456:G:C"), data.frame(CHROM = "chrX", POS = 456, REF = "G", ALT = "C", stringsAsFactors = FALSE))

  # Special characters for REF and ALT
  expect_equal(read_mut("4:111:_:-"), data.frame(CHROM = "4", POS = 111, REF = "_", ALT = "-", stringsAsFactors = FALSE))

  # NA in REF of ALT
  expect_equal(read_mut("2:678:NA:T"), data.frame(CHROM = "2", POS = 678, REF = "NA", ALT = "T", stringsAsFactors = FALSE))

  # Multiallelic cases
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

  # Invalid case to see error
  # Missing Ref or Alt
  expect_error(
    read_mut("chr1:123:A"),
    "The parameter 'mut' ('chr1:123:A') is not in the expected format",
    fixed = TRUE
  )
  # Add extra parameter
  expect_error(
    read_mut("chr1:123:A:T:extra"),
    "The parameter 'mut' ('chr1:123:A:T:extra') is not in the expected format",
    fixed = TRUE
  )

  # Invalid_format
  expect_error(
    read_mut("invalid_format"),
    "The parameter 'mut' ('invalid_format') is not in the expected format",
    fixed = TRUE
  )
})

test_that("read_mut works with VCF files", {
  # ----------------------------------------------------------------------------------------------
  # SETUP: Synthetic VCF Content matching
  # ----------------------------------------------------------------------------------------------
  vcf_content <- c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t12345\t.\tA\tT\t.\t.\t.",
    "2\t67890\t.\tGT\tC\t.\t.\t.",
    "chr3\t101112\t.\tA.\tA-\t.\t.\t.",
    "chr22\t56789\t.\tNA\tNA\t.\t.\t.",
    "22\t56789\t.\t.\t-\t.\t.\t.",
    "chrX\t54321\t.\tA\t.\t.\t.\t.",
    "chrY\t99999\t.\t.\tC\t.\t.\t.",
    "chr4\t1234\t.\t-\tA\t.\t.\t.",
    "chr5\t5678\t.\t.\t.\t.\t.\t.",
    "chr6\t91011\t.\tG\tNA\t.\t.\t.",
    "chr8\t151617\t.\tG\tA,T\t.\t.\t.",
    "chrY\tNA\t.\tT\tC\t.\t.\t.",
    "42453\t181920\t.\tAC\t-\t.\t.\t.",
    ".\t181920\t.\tAC\t-\t.\t.\t."
  )

  vcf_file <- tempfile(fileext = ".vcf")
  writeLines(vcf_content, vcf_file)
  on.exit(unlink(vcf_file))

  # Case 1: Standard VCF
  captured_warnings <- testthat::capture_warnings(
    df_vcf_test <- read_mut(vcf_file)
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

  # Check warnings
  expect_true(any(grepl("non-integer values in the POS column", captured_warnings)))

  # Check results
  # We sort to ensure stable comparison if row order varies
  expect_equal(df_vcf_test, expected_results)

  # Case 2: Compressed VCF
  vcf_file_gz <- tempfile(fileext = ".vcf.gz")
  # Use gzfile to write compressed
  gz_con <- gzfile(vcf_file_gz, "wt")
  writeLines(vcf_content, gz_con)
  close(gz_con)
  on.exit(unlink(vcf_file_gz), add = TRUE)

  captured_warnings_gz <- testthat::capture_warnings(
    df_vcf_test2 <- read_mut(vcf_file_gz)
  )

  expect_true(any(grepl("non-integer values in the POS column", captured_warnings_gz)))
  expect_equal(df_vcf_test2, expected_results)
})

test_that("read_mut works with TSV files", {
  # ----------------------------------------------------------------------------------------------
  # SETUP: Synthetic TSV Content
  # ----------------------------------------------------------------------------------------------
  # TSV Header: CHROM POS REF ALT
  tsv_content <- c(
    "CHROM\tPOS\tREF\tALT",
    "chr1\t12345\tA\tT",
    "2\t67890\tGT\tC",
    "chr3\t101112\tA.\tA-",
    "chr22\t56789\tNA\tNA",
    "22\t56789\t.\t-",
    "chrX\t54321\tA\t-",
    "chrY\t99999\t\tC",
    "chr4\t1234\t-\tA",
    "chr5\t5678\t\t",
    "chr6\t91011\tG\tNA",
    "chr8\t151617\tG\tA,T",
    "chrX\tNA\tA\t-",
    "chrY\tNA\tT\tC",
    "42453\t181920\tAC\t-",
    "\t181920\tAC\t-"
  )

  tsv_file <- tempfile(fileext = ".tsv")
  writeLines(tsv_content, tsv_file)
  on.exit(unlink(tsv_file))

  # Case 1: Standard TSV
  captured_warnings <- testthat::capture_warnings(
    df_tsv_test <- read_mut(tsv_file)
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

  # Check warnings
  expect_true(any(grepl("non-integer values in the POS column", captured_warnings)))

  # Check results
  expect_equal(df_tsv_test, expected_results3)

  # Case 2: Compressed TSV
  tsv_file_gz <- tempfile(fileext = ".tsv.gz")
  gz_con <- gzfile(tsv_file_gz, "wt")
  writeLines(tsv_content, gz_con)
  close(gz_con)
  on.exit(unlink(tsv_file_gz), add = TRUE)

  captured_warnings_gz <- testthat::capture_warnings(
    df_tsv_test2 <- read_mut(tsv_file_gz)
  )

  expect_true(any(grepl("non-integer values in the POS column", captured_warnings_gz)))
  expect_equal(df_tsv_test2, expected_results3)
})

test_that("read_mut detects multiallelic REF in VCF/TSV", {
  # ----------------------------------------------------------------------------------------------
  # SETUP: VCF with Multiallelic REF (invalid)
  # ----------------------------------------------------------------------------------------------
  vcf_content <- c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA,G\tT\t.\t.\t."
  )
  vcf_file <- tempfile(fileext = ".vcf")
  writeLines(vcf_content, vcf_file)
  on.exit(unlink(vcf_file))

  captured_warnings <- testthat::capture_warnings(
    expect_error(read_mut(vcf_file), "After reading and filtering, the mutation information is empty")
  )
  expect_true(any(grepl("REF can not be multiallelic", captured_warnings)))
})
