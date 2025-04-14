test_that("Read bam", {
  bam_test <- system.file("testdata/bam/", "chr17_examples.sorted.bam", package = "fRagmentomics")

  # Call the read bam function
  result <- read_bam(
    bam         = bam_test,
    chr         = "chr17",
    pos         = 7578063,
    neg_offset  = 20,
    pos_offset  = 20,
    flag_keep   = 0x03,
    flag_remove = 0x900
  )

  # Create the fake dataframe to compare
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
    RNAME = c("chr17", "chr17", "chr17", "chr17", "chr17", "chr17"),
    POS = c(7578027, 7578047, 7578058, 7578063, 7578063, 7578078),
    TLEN = c(196, 160, 150, -160, -150, -196),
    MAPQ = rep(60, 6),
    CIGAR = c("144M", "144M", "20M1D124M", "144M", "15M1D129M", "119M1D25M"),
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
    )
  )

  # Compare it's the expected result
  expect_equal(result, expected_df)
})

test_that("Read bam", {
  bam_test <- system.file("testdata/bam/", "cfdna-test-01_chr17_7576000_7579000.bam", package = "fRagmentomics")

  result <- read_bam(
    bam         = bam_test,
    chr         = "chr17",
    pos         = 2191,
    neg_offset  = -1000,
    pos_offset  = 1000,
    flag_keep   = 0x03,
    flag_remove = 0x900
  )

  expected_cols <- c(
    "QNAME", "FLAG", "RNAME", "POS",
    "TLEN", "MAPQ", "CIGAR", "SEQ",
    "QUAL"
  )

  expect_true(
    all(expected_cols %in% colnames(result))
  )
})

test_that("Read bam empty", {
  bam_test <- system.file("testdata/bam/", "cfdna-test-01_chr17_7576000_7579000.bam", package = "fRagmentomics")

  expect_error(
    read_bam(
      bam         = bam_test,
      chr         = "chr17",
      pos         = 15168135,
      neg_offset  = -1000,
      pos_offset  = 1000,
      flag_keep   = 0x03,
      flag_remove = 0x900
    ),
    "The final BAM data frame is empty. No reads match the criteria."
  )
})

test_that("Read bam missing columns", {
  bam_test <- system.file("testdata/bam/", "chr17_examples.sorted.bam", package = "fRagmentomics")

  result <- read_bam(
    bam         = bam_test,
    chr         = "chr17",
    pos         = 7578063,
    neg_offset  = -20,
    pos_offset  = 20,
    flag_keep   = 0x03,
    flag_remove = 0x900
  )

  # Simulation when deletion of a column
  result <- result[, -which(names(result) == "QUAL")]

  expect_error(
    {
      required_columns <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "SEQ", "QUAL")
      missing_columns <- setdiff(required_columns, colnames(result))
      if (length(missing_columns) > 0) {
        stop(paste("The final BAM data frame is missing the following columns:", paste(missing_columns, collapse = ", ")))
      }
    },
    regexp = "missing the following columns: QUAL"
  )
})
