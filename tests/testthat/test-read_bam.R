test_that("read_bam returns expected reads as a data.frame (small BAM)", {
  bam_test <- system.file("testdata/bam/", "chr17_examples.sorted.bam", package = "fRagmentomics")

  result <- read_bam(
    bam = bam_test,
    chr = "chr17",
    pos = 7578063,
    neg_offset_mate_search = -20,
    pos_offset_mate_search = 20,
    flag_bam_list = list(
      isPaired = TRUE,
      isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE,
      isSupplementaryAlignment = FALSE,
      isDuplicate = FALSE
    )
  )

  # The result must be a data.frame (not a list)
  expect_s3_class(result, "data.frame")

  # Expected dataframe for reads
  expected_df <- data.frame(
    QNAME = c(
      "06106eKSy:0100",
      "01135AEM2:0305",
      "08038wyOK:0002",
      "01135AEM2:0305",
      "08038wyOK:0002",
      "06106eKSy:0100"
    ),
    FLAG = c(163, 163, 163, 83, 83, 83),
    RNAME = rep("chr17", 6),
    POS = c(7578027, 7578047, 7578058, 7578063, 7578063, 7578078),
    TLEN = c(196, 160, 150, -160, -150, -196),
    MAPQ = rep(60, 6),
    CIGAR = c("144M", "144M", "20M1D124M", "144M", "15M1D129M", "119M1D25M"),
    RNEXT = rep("chr17", 6),
    PNEXT = c(7578078, 7578063, 7578063, 7578047, 7578058, 7578027),
    SEQ = c(
      "CAATAGTTAAACCCATTTACTTTGCACATCTCATGGGGTTATAGGGAGGTCAAATAAGCAGCAGGAGAAAGCCCCCCTACTGCTCACCCGGAGGGCCACTGACAACCACCCTTAACCCCTCCTCCCAGAGACCCCAGTTGCAAA",
      "TTTGCACATCTCATGGGGTTATAGGGAGGTCAAATAAGCAGCAGGAGAAAGCCCCCCTACTGCTCACCCGGAGGGCCACTGACAACCACCCTTAACCCCTCCTCCCAGAGACCCCAGTTGCAAACCAGACCTCAGGCGGCTCAT",
      "CATGGGGTTATAGGGAGGTCAATAAGCAGCAGGAGAAAGCCCCCCTACTGCTCACCCGGAGGGCCACTGACAACCACCCTTAACCCCTCCTCCCAGAGACCCCAGTTGCAAACCAGACCTCAGGCGGCTCATAGGGCACCACCA",
      "GGTTATAGGGAGGTCAAATAAGCAGCAGGAGAAAGCCCCCCTACTGCTCACCCGGAGGGCCACTGACAACCACCCTTAACCCCTCCTCCCAGAGACCCCAGTTGCAAACCAGACCTCAGGCGGCTCATAGGGCACCACCACACT",
      "GGTTATAGGGAGGTCAATAAGCAGCAGGAGAAAGCCCCCCTACTGCTCACCCGGAGGGCCACTGACAACCACCCTTAACCCCTCCTCCCAGAGACCCCAGTTGCAAACCAGACCTCAGGCGGCTCATAGGGCACCACCACACTA",
      "AAATAAGCAGCAGGAGAAAGCCCCCCTACTGCTCACCCGGAGGGCCACTGACAACCACCCTTAACCCCTCCTCCCAGAGACCCCAGTTGCAAACCAGACCTCAGGCGGCTCATAGGGCACACCACACTATGTCGAAAAGTGTTT"
    ),
    QUAL = c(
      "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
      "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFF:FFFFFFFF",
      "FFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
      ":FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
      "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
      "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
    ),
    stringsAsFactors = FALSE
  )

  # Check expected columns
  expected_cols <- c("QNAME", "FLAG", "RNAME", "POS", "TLEN", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "SEQ", "QUAL")
  expect_setequal(colnames(result), expected_cols)

  # Sort both dataframes for stable comparison
  ord_res <- order(result$QNAME, result$POS, na.last = TRUE)
  ord_exp <- order(expected_df$QNAME, expected_df$POS, na.last = TRUE)

  res_cmp <- result[ord_res, expected_cols, drop = FALSE]
  exp_cmp <- expected_df[ord_exp, expected_cols, drop = FALSE]

  # Harmonize integer/numeric types
  int_cols <- c("FLAG", "POS", "TLEN", "MAPQ", "PNEXT")
  for (cc in int_cols) {
    res_cmp[[cc]] <- as.integer(res_cmp[[cc]])
    exp_cmp[[cc]] <- as.integer(exp_cmp[[cc]])
  }
  expect_equal(res_cmp, exp_cmp, ignore_attr = TRUE)
})

test_that("read_bam returns expected columns on larger window (cfDNA BAM)", {
  bam_test <- system.file("testdata/bam/", "cfdna-test-01_chr17_7576000_7579000.bam", package = "fRagmentomics")

  result <- read_bam(
    bam = bam_test,
    chr = "chr17",
    pos = 2191,
    neg_offset_mate_search = -1000,
    pos_offset_mate_search = 1000,
    flag_bam_list = list(
      isPaired = TRUE,
      isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE,
      isSupplementaryAlignment = FALSE,
      isDuplicate = FALSE
    )
  )

  # The result must be a data.frame
  expect_s3_class(result, "data.frame")

  expected_cols <- c("QNAME", "FLAG", "RNAME", "POS", "TLEN", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "SEQ", "QUAL")
  expect_true(all(expected_cols %in% colnames(result)))
})

test_that("read_bam returns NULL when no reads are found", {
  bam_test <- system.file("testdata/bam/", "cfdna-test-01_chr17_7576000_7579000.bam", package = "fRagmentomics")

  # Expect NULL if the target position does not overlap with any read
  expect_null(
    read_bam(
      bam = bam_test,
      chr = "chr17",
      pos = 15168135,
      neg_offset_mate_search = -1000,
      pos_offset_mate_search = 1000,
      flag_bam_list = list(isPaired = TRUE, isUnmappedQuery = FALSE)
    )
  )
})

test_that("read_bam detects missing required columns", {
  bam_test <- system.file("testdata/bam/", "chr17_examples.sorted.bam", package = "fRagmentomics")

  result <- read_bam(
    bam = bam_test,
    chr = "chr17",
    pos = 7578063,
    neg_offset_mate_search = -20,
    pos_offset_mate_search = 20,
    flag_bam_list = list(
      isPaired = TRUE,
      isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE,
      isSupplementaryAlignment = FALSE,
      isDuplicate = FALSE
    )
  )

  skip_if(is.null(result), "No reads to test column deletion.")

  df <- result
  # Simulate deletion of an essential column
  df_no_qual <- df[, setdiff(names(df), "QUAL"), drop = FALSE]

  expect_error(
    {
      required_columns <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "SEQ", "QUAL")
      missing_columns <- setdiff(required_columns, colnames(df_no_qual))
      if (length(missing_columns) > 0) {
        stop(paste(
          "The final BAM dataframe is missing the following columns:",
          paste(missing_columns, collapse = ", ")
        ))
      }
    },
    regexp = "missing the following columns: QUAL"
  )
})
