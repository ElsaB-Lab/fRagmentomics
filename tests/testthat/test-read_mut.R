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
  # SETUP: Valid VCF content (SNV, deletion, multi-allelic)
  # ----------------------------------------------------------------------------------------------
  vcf_content <- c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t12345\t.\tA\tT\t.\t.\t.",
    "2\t67890\t.\tGAT\tG\t.\t.\t.",
    "chr8\t151617\t.\tG\tA,T\t.\t.\t."
  )

  vcf_file <- tempfile(fileext = ".vcf")
  writeLines(vcf_content, vcf_file)
  on.exit(unlink(vcf_file))

  expected_results <- data.frame(
    CHROM = c("chr1", "2", "chr8", "chr8"),
    POS   = c(12345L, 67890L, 151617L, 151617L),
    REF   = c("A", "GAT", "G", "G"),
    ALT   = c("T", "G", "A", "T"),
    stringsAsFactors = FALSE
  )

  # Case 1: Plain VCF
  df_vcf_test <- read_mut(vcf_file)
  expect_equal(df_vcf_test, expected_results)

  # Case 2: Compressed VCF (bgzip + tabix index required by VariantAnnotation)
  vcf_gz <- tempfile(fileext = ".vcf.gz")
  Rsamtools::bgzip(vcf_file, dest = vcf_gz, overwrite = TRUE)
  Rsamtools::indexTabix(vcf_gz, format = "vcf")
  on.exit(unlink(c(vcf_gz, paste0(vcf_gz, ".tbi"))), add = TRUE)

  df_vcf_gz <- read_mut(vcf_gz)
  expect_equal(df_vcf_gz, expected_results)
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

test_that("read_mut detects multiallelic REF in TSV", {
  # ----------------------------------------------------------------------------------------------
  # SETUP: TSV with Multiallelic REF (invalid)
  # Note: VCF spec forbids multiallelic REF, so this is tested via TSV input
  # ----------------------------------------------------------------------------------------------
  tsv_content <- c(
    "CHROM\tPOS\tREF\tALT",
    "chr1\t100\tA,G\tT"
  )
  tsv_file <- tempfile(fileext = ".tsv")
  writeLines(tsv_content, tsv_file)
  on.exit(unlink(tsv_file))

  captured_warnings <- testthat::capture_warnings(
    expect_error(read_mut(tsv_file), "After reading and filtering, the mutation information is empty")
  )
  expect_true(any(grepl("REF can not be multiallelic", captured_warnings)))
})

test_that("parser_vcf handles header-only VCF and corrupted compressed files", {
  # ----------------------------------------------------------------------------------------------
  # SETUP: Edge cases for parser_vcf low-level error branches
  # ----------------------------------------------------------------------------------------------

  # Header-only VCF (no variant rows) -> length(vcf) == 0 branch
  vcf_header_only <- c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  )
  vcf_empty <- tempfile(fileext = ".vcf")
  writeLines(vcf_header_only, vcf_empty)
  on.exit(unlink(vcf_empty))

  expect_error(
    read_mut(vcf_empty),
    "The VCF file does not contain any mutation data."
  )

  # Corrupted .vcf.gz file -> tryCatch error branch in parser_vcf
  vcf_gz_bad <- tempfile(fileext = ".vcf.gz")
  writeBin(charToRaw("this is not valid gzip content"), vcf_gz_bad)
  on.exit(unlink(vcf_gz_bad), add = TRUE)

  expect_error(
    read_mut(vcf_gz_bad),
    "Failed to read VCF file"
  )
})

test_that("parser_tsv handles missing columns, unreadable files, and all-empty rows", {
  # ----------------------------------------------------------------------------------------------
  # SETUP: Edge cases for parser_tsv low-level error branches
  # ----------------------------------------------------------------------------------------------

  # TSV with missing required columns -> ncol(df_tsv) != 4 branch
  # cols_only() returns only matched columns, so missing headers reduce ncol
  tsv_missing_cols <- c(
    "CHROM\tPOS\tREF",
    "chr1\t123\tA"
  )
  tsv_file_missing <- tempfile(fileext = ".tsv")
  writeLines(tsv_missing_cols, tsv_file_missing)
  on.exit(unlink(tsv_file_missing))

  captured_warnings_missing <- testthat::capture_warnings(
    expect_error(
      read_mut(tsv_file_missing),
      "does not match the expected columns"
    )
  )
  expect_true(any(grepl("don't match the column names", captured_warnings_missing)))

  # Non-existent TSV file -> tryCatch error branch in parser_tsv
  expect_error(
    parser_tsv(tempfile(fileext = ".tsv")),
    "Failed to read the TSV file"
  )

  # TSV with only fully empty rows -> nrow(df_tsv) == 0 after filtering branch
  tsv_all_empty <- c(
    "CHROM\tPOS\tREF\tALT",
    "\t\t\t",
    "\t\t\t"
  )
  tsv_file_empty <- tempfile(fileext = ".tsv")
  writeLines(tsv_all_empty, tsv_file_empty)
  on.exit(unlink(tsv_file_empty), add = TRUE)

  expect_error(
    read_mut(tsv_file_empty),
    "The TSV file does not contain any mutation data."
  )
})

test_that("parser_str handles trailing colon, non-integer POS, and length != 4", {
  # ----------------------------------------------------------------------------------------------
  # SETUP: Edge cases for parser_str low-level branches
  # Note: the trailing-colon and length != 4 branches are unreachable via read_mut because
  # ----------------------------------------------------------------------------------------------

  # Trailing colon -> mut is patched to append "-" before splitting
  result_trailing <- parser_str("chr1:123:A:")
  expect_equal(
    result_trailing,
    data.frame(CHROM = "chr1", POS = 123L, REF = "A", ALT = "-", stringsAsFactors = FALSE)
  )

  # length(parts) != 4 after split -> error branch in parser_str
  expect_error(
    parser_str("chr1:123:A:T:extra"),
    "not in the expected format"
  )

  # Non-integer POS in string format -> warning + NA
  expect_warning(
    result_na_pos <- read_mut("chr1:abc:A:T"),
    "not a valid integer"
  )
  expect_true(is.na(result_na_pos$POS))
  expect_equal(result_na_pos$CHROM, "chr1")
  expect_equal(result_na_pos$REF,   "A")
  expect_equal(result_na_pos$ALT,   "T")
})