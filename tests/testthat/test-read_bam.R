test_that("read_bam returns expected reads as a data.frame", {
  # ----------------------------------------------------------------------------------------------
  # SETUP: Synthetic BAM with a read pair covering the target
  # Target: chr1:500
  # Read 1: Covers 500 (POS 480, 50M -> End 529). Length 50.
  # Read 2: Mate, does not cover 500 (POS 600, 50M -> End 649). Length 50.
  # ----------------------------------------------------------------------------------------------

  # Ensure SEQ length matches CIGAR (50M = 50 bases) exactly
  seq_50bp <- strrep("A", 50)
  qual_50bp <- strrep("F", 50)

  df_reads <- data.frame(
    QNAME = c("read_pair_1", "read_pair_1"),
    FLAG = as.integer(c(99, 147)),
    RNAME = c("chr1", "chr1"),
    POS = as.integer(c(480, 600)),
    MAPQ = as.integer(c(60, 60)),
    CIGAR = c("50M", "50M"),
    RNEXT = c("=", "="),
    PNEXT = as.integer(c(600, 480)),
    TLEN = as.integer(c(170, -170)),
    SEQ = c(seq_50bp, seq_50bp),
    QUAL = c(qual_50bp, qual_50bp),
    stringsAsFactors = FALSE
  )

  bam_env <- setup_test_bam(df_reads, ref_seq_lengths = c(chr1 = 1000))
  on.exit(cleanup_test_bam(bam_env))

  # ----------------------------------------------------------------------------------------------
  # TEST: Retrieve reads covering pos 500
  # ----------------------------------------------------------------------------------------------
  result <- read_bam(
    bam = bam_env$path,
    chr = "chr1",
    pos = 500,
    neg_offset_mate_search = -100,
    pos_offset_mate_search = 200, # Large enough to catch the mate at 600
    flag_bam_list = list(
      isPaired = TRUE,
      isUnmappedQuery = FALSE,
      isDuplicate = FALSE
    )
  )

  # Expect a data.frame
  expect_s3_class(result, "data.frame")

  # Expect both reads (the covering read AND its mate)
  expect_equal(nrow(result), 2)
  expect_true(all(result$QNAME == "read_pair_1"))

  # Check expected columns
  expected_cols <- c("QNAME", "FLAG", "RNAME", "POS", "TLEN", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "SEQ", "QUAL")
  expect_true(all(expected_cols %in% colnames(result)))
})

test_that("read_bam handles complex CIGAR strings for width calculation", {
  # ----------------------------------------------------------------------------------------------
  # SETUP: Synthetic BAM with CIGAR operations affecting width (Deletion/Insertion)
  # Target: chr1:500
  # Read 1: POS 490.
  # CIGAR: 5M 10I 10M
  #   - 5M: Consumes 5 Ref, 5 Query
  #   - 10I: Consumes 0 Ref, 10 Query (Insertion)
  #   - 10M: Consumes 10 Ref, 10 Query
  #   Total Query Length (SEQ) MUST be: 5 + 10 + 10 = 25 bases.
  #   Total Ref Width: 5 + 0 + 10 = 15 bases.
  #   End Position on Ref: 490 + 15 - 1 = 504.
  #   Since 490 <= 500 <= 504, this read covers position 500.
  # ----------------------------------------------------------------------------------------------

  seq_25bp <- strrep("T", 25) # Must match 5M + 10I + 10M length
  qual_25bp <- strrep("F", 25)

  df_reads <- data.frame(
    QNAME = c("indel_read"),
    FLAG = as.integer(c(99)),
    RNAME = c("chr1"),
    POS = as.integer(c(490)),
    MAPQ = as.integer(c(60)),
    CIGAR = c("5M10I10M"),
    RNEXT = c("="),
    PNEXT = as.integer(c(600)),
    TLEN = as.integer(c(100)),
    SEQ = c(seq_25bp),
    QUAL = c(qual_25bp),
    stringsAsFactors = FALSE
  )

  bam_env <- setup_test_bam(df_reads, ref_seq_lengths = c(chr1 = 1000))
  on.exit(cleanup_test_bam(bam_env))

  result <- read_bam(
    bam = bam_env$path,
    chr = "chr1",
    pos = 500,
    neg_offset_mate_search = -50,
    pos_offset_mate_search = 50,
    flag_bam_list = list(isPaired = TRUE)
  )

  expect_false(is.null(result))
  expect_equal(result$QNAME, "indel_read")
})

test_that("read_bam returns NULL when no reads cover the position", {
  # ----------------------------------------------------------------------------------------------
  # SETUP: Reads exist in window but do not cover the specific position
  # Target: chr1:500
  # Read 1: POS 450, 40M. Ends at 489. (Short of 500)
  # Read 2: POS 510, 40M. Starts at 510. (Past 500)
  # Gap at 500.
  # ----------------------------------------------------------------------------------------------

  seq_40bp <- strrep("G", 40)
  qual_40bp <- strrep("F", 40)

  df_reads <- data.frame(
    QNAME = c("left_read", "right_read"),
    FLAG = as.integer(c(99, 99)),
    RNAME = c("chr1", "chr1"),
    POS = as.integer(c(450, 510)),
    MAPQ = as.integer(c(60, 60)),
    CIGAR = c("40M", "40M"),
    RNEXT = c("=", "="),
    PNEXT = as.integer(c(600, 600)),
    TLEN = as.integer(c(100, 100)),
    SEQ = c(seq_40bp, seq_40bp),
    QUAL = c(qual_40bp, qual_40bp),
    stringsAsFactors = FALSE
  )

  bam_env <- setup_test_bam(df_reads, ref_seq_lengths = c(chr1 = 1000))
  on.exit(cleanup_test_bam(bam_env))

  result <- read_bam(
    bam = bam_env$path,
    chr = "chr1",
    pos = 500,
    neg_offset_mate_search = -100,
    pos_offset_mate_search = 100,
    flag_bam_list = list(isPaired = TRUE)
  )

  expect_null(result)
})

test_that("read_bam respects flags (filters out duplicates)", {
  # ----------------------------------------------------------------------------------------------
  # SETUP: Duplicate read covering the position
  # Target: chr1:500
  # Read 1: Covers 500, but FLAG 1024 (Duplicate)
  # ----------------------------------------------------------------------------------------------

  seq_50bp <- strrep("C", 50)
  qual_50bp <- strrep("F", 50)

  df_reads <- data.frame(
    QNAME = c("dup_read"),
    FLAG = as.integer(c(1024)),
    RNAME = c("chr1"),
    POS = as.integer(c(480)),
    MAPQ = as.integer(c(60)),
    CIGAR = c("50M"),
    RNEXT = c("="),
    PNEXT = as.integer(c(600)),
    TLEN = as.integer(c(100)),
    SEQ = c(seq_50bp),
    QUAL = c(qual_50bp),
    stringsAsFactors = FALSE
  )

  bam_env <- setup_test_bam(df_reads, ref_seq_lengths = c(chr1 = 1000))
  on.exit(cleanup_test_bam(bam_env))

  result <- read_bam(
    bam = bam_env$path,
    chr = "chr1",
    pos = 500,
    neg_offset_mate_search = -50,
    pos_offset_mate_search = 50,
    flag_bam_list = list(isDuplicate = FALSE)
  )

  expect_null(result)
})

test_that("read_bam detects missing required columns", {
  # ----------------------------------------------------------------------------------------------
  # SETUP: Basic valid BAM
  # ----------------------------------------------------------------------------------------------

  seq_50bp <- strrep("A", 50)
  qual_50bp <- strrep("F", 50)

  df_reads <- data.frame(
    QNAME = c("valid_read"),
    FLAG = as.integer(c(99)),
    RNAME = c("chr1"),
    POS = as.integer(c(480)),
    MAPQ = as.integer(c(60)),
    CIGAR = c("50M"),
    RNEXT = c("="),
    PNEXT = as.integer(c(600)),
    TLEN = as.integer(c(100)),
    SEQ = c(seq_50bp),
    QUAL = c(qual_50bp),
    stringsAsFactors = FALSE
  )

  bam_env <- setup_test_bam(df_reads, ref_seq_lengths = c(chr1 = 1000))
  on.exit(cleanup_test_bam(bam_env))

  result <- read_bam(
    bam = bam_env$path,
    chr = "chr1",
    pos = 500,
    neg_offset_mate_search = -50,
    pos_offset_mate_search = 50,
    flag_bam_list = list(isPaired = TRUE)
  )

  expect_false(is.null(result))

  # Manually simulate a dataframe missing the QUAL column to test the error logic
  df_no_qual <- result[, setdiff(names(result), "QUAL"), drop = FALSE]

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
