test_that("Test function remove softclip", {
  # Test 1: Correctly removes soft clips from both ends
  read_stats_input <- list(
    CIGAR = "10S80M10S",
    SEQ = paste0(
      paste(rep("A", 10), collapse = ""),
      paste(rep("G", 80), collapse = ""),
      paste(rep("C", 10), collapse = "")
    ),
    QUAL = paste(rep("F", 100), collapse = "")
  )
  result <- remove_softclip(read_stats_input)

  # Assert: Define expected output and check if the result matches
  expected_output <- list(
    SEQ = paste(rep("G", 80), collapse = ""),
    QUAL = substr(read_stats_input$QUAL, 11, 90),
    CIGAR = "80M",
    read_length = 80
  )

  expect_equal(result, expected_output)

  # Test 2: Correctly removes 5p soft clips only
  read_stats_input <- list(
    CIGAR = "5S95M",
    SEQ = paste0(
      paste(rep("T", 5), collapse = ""),
      paste(rep("G", 95), collapse = "")
    ),
    QUAL = paste(rep("F", 100), collapse = "")
  )
  result <- remove_softclip(read_stats_input)

  # Assert
  expected_output <- list(
    SEQ = paste(rep("G", 95), collapse = ""),
    QUAL = substr(read_stats_input$QUAL, 6, 100),
    CIGAR = "95M",
    read_length = 95
  )

  expect_equal(result, expected_output)

  # Test 3: Correctly removes 3p soft clips only
  read_stats_input <- list(
    CIGAR = "92M8S",
    SEQ = paste0(
      paste(rep("G", 92), collapse = ""),
      paste(rep("C", 8), collapse = "")
    ),
    QUAL = paste(rep("F", 100), collapse = "")
  )
  result <- remove_softclip(read_stats_input)

  # Assert
  expected_output <- list(
    SEQ = paste(rep("G", 92), collapse = ""),
    QUAL = substr(read_stats_input$QUAL, 1, 92),
    CIGAR = "92M",
    read_length = 92
  )

  expect_equal(result, expected_output)

  # Test 4: Returns the original data when there are no soft clips
  read_stats_input <- list(
    SEQ   = paste(rep("G", 100), collapse = ""),
    QUAL  = paste(rep("F", 100), collapse = ""),
    CIGAR = "100M"
  )
  result <- remove_softclip(read_stats_input)

  # Assert: The result should be identical to the input
  expected_output <- list(
    SEQ = paste(rep("G", 100), collapse = ""),
    QUAL = paste(rep("F", 100), collapse = ""),
    CIGAR = "100M",
    read_length = 100
  )
  expect_equal(result, expected_output)


  # Test 5: Returns empty strings when the entire read is soft-clipped
  read_stats_input <- list(
    CIGAR = "150S",
    SEQ   = paste(rep("N", 150), collapse = ""),
    QUAL  = paste(rep("#", 150), collapse = "")
  )
  result <- remove_softclip(read_stats_input)

  # Assert: All fields should be empty strings
  expected_output <- list(
    SEQ = "",
    QUAL = "",
    CIGAR = "",
    read_length = 0
  )

  expect_equal(result, expected_output)


  # Test 6: Correctly handles complex CIGARs with insertions/deletions
  read_stats_input <- list(
    CIGAR = "10S40M5I45M5S",
    SEQ   = paste(rep("A", 100), collapse = ""),
    QUAL  = paste(rep("F", 100), collapse = "")
  )
  result <- remove_softclip(read_stats_input)

  # Assert
  expected_output <- list(
    SEQ = substr(read_stats_input$SEQ, 11, 95), # 100 - 10 - 5 = 85 chars long
    QUAL = substr(read_stats_input$QUAL, 11, 95),
    CIGAR = "40M5I45M",
    read_length = 85
  )

  expect_equal(result, expected_output)

  # Test 7: Correctly ignores hard clips and only removes soft clips
  read_stats_input <- list(
    CIGAR = "5H10S80M5S5H",
    SEQ   = paste(rep("A", 95), collapse = ""),
    QUAL  = paste(rep("F", 95), collapse = "")
  )
  result <- remove_softclip(read_stats_input)

  # Assert: The 'S' parts should be gone, but 'H' should remain
  expected_output <- list(
    SEQ = substr(read_stats_input$SEQ, 11, 90), # 95 - 10 - 5 = 80 chars
    QUAL = substr(read_stats_input$QUAL, 11, 90),
    CIGAR = "5H80M5H",
    read_length = 80
  )

  expect_equal(result, expected_output)
})

test_that("remove_softclip trims both ends and yields the expected fragment size", {
  # Synthetic paired-end fragment with soft clipping:
  #  - 5' read (forward):    5S20M starting at POS=100  -> CIGAR '5S20M'
  #  - 3' read (reverse):    20M5S whose aligned 20M spans 181..200 -> CIGAR '20M5S'
  #    (leftmost POS for the reverse read is 181 so that the aligned end is at 200)
  #
  # Flags: 99 (R1 forward), 147 (R2 reverse) are standard proper-pair flags.
  df_sam <- data.frame(
    QNAME = c("frag_soft", "frag_soft"),
    FLAG = c(99L, 147L),
    RNAME = c("chr1", "chr1"),
    POS = c(100L, 181L),
    MAPQ = c(60L, 60L),
    CIGAR = c("5S20M", "20M5S"),
    RNEXT = c("=", "="),
    PNEXT = c(181L, 100L),
    TLEN = c(101L, -101L),
    SEQ = c(
      paste0(strrep("N", 5), strrep("A", 20)),
      paste0(strrep("T", 20), strrep("N", 5))
    ),
    QUAL = c(strrep("I", 25), strrep("I", 25)),
    stringsAsFactors = FALSE
  )

  # Minimal locus/context to satisfy downstream calls; we don't assert on mutation status.
  chr <- "chr1"
  pos <- 150L
  ref <- "A"
  alt <- "T"
  fasta_seq <- list(chr = "chr1", start = 1L, end = 1000L, seq = paste(rep("A", 1000), collapse = ""))

  # --- Run WITHOUT trimming ----------------------------------------------------
  res_no_trim <- extract_fragment_features(
    df_sam = df_sam,
    fragment_name = "frag_soft",
    sample_id = "sampleX",
    chr = chr, pos = pos, ref = ref, alt = alt,
    report_softclip = TRUE,
    report_5p_3p_bases_fragment = 0L, # avoid extra work; not needed here
    remove_softclip = FALSE,
    report_bam_info = TRUE,
    fasta_seq = fasta_seq,
    input_mutation_info = "chr1:150:A>T"
  )

  # Expect original CIGARs preserved
  expect_identical(res_no_trim$CIGAR_5p, "5S20M")
  expect_identical(res_no_trim$CIGAR_3p, "20M5S")

  # Soft-clip counts (computed on the returned CIGARs)
  expect_equal(
    c(res_no_trim$Nb_Fragment_Bases_Softclip_5p, res_no_trim$Nb_Fragment_Bases_Softclip_3p),
    c(5L, 5L)
  )

  # Also compute the "no trim" fragment size using the package helper for comparison
  rs5_orig <- list(
    FLAG = df_sam$FLAG[1], MAPQ = df_sam$MAPQ[1], TLEN = df_sam$TLEN[1],
    CIGAR = df_sam$CIGAR[1], POS = df_sam$POS[1],
    SEQ = df_sam$SEQ[1], QUAL = df_sam$QUAL[1], read_length = nchar(df_sam$SEQ[1])
  )
  rs3_orig <- list(
    FLAG = df_sam$FLAG[2], MAPQ = df_sam$MAPQ[2], TLEN = df_sam$TLEN[2],
    CIGAR = df_sam$CIGAR[2], POS = df_sam$POS[2],
    SEQ = df_sam$SEQ[2], QUAL = df_sam$QUAL[2], read_length = nchar(df_sam$SEQ[2])
  )
  size_no_trim_expected <- fRagmentomics:::get_fragment_size(rs5_orig, rs3_orig)
  expect_identical(res_no_trim$Fragment_Size, as.integer(size_no_trim_expected))

  # --- Run WITH trimming -------------------------------------------------------
  res_trim <- extract_fragment_features(
    df_sam = df_sam,
    fragment_name = "frag_soft",
    sample_id = "sampleX",
    chr = chr, pos = pos, ref = ref, alt = alt,
    report_bam_info = TRUE,
    report_softclip = TRUE,
    report_5p_3p_bases_fragment = 0L,
    remove_softclip = TRUE,
    fasta_seq = fasta_seq,
    input_mutation_info = "chr1:150:A>T"
  )

  # We keep the input information for BAM information
  expect_identical(res_trim$CIGAR_5p, "5S20M")
  expect_identical(res_trim$CIGAR_3p, "20M5S")

  # Soft-clip counts must be zero on both fragment ends after trimming
  # This feature is calculated after the softclip
  expect_equal(
    c(res_trim$Nb_Fragment_Bases_Softclip_5p, res_trim$Nb_Fragment_Bases_Softclip_3p),
    c(0L, 0L)
  )

  # The fragment size after trimming must match what the helpers compute on
  # the trimmed reads (and cannot be larger than the untrimmed size).
  rs5_trim <- utils::modifyList(rs5_orig, fRagmentomics:::remove_softclip(rs5_orig))
  rs3_trim <- utils::modifyList(rs3_orig, fRagmentomics:::remove_softclip(rs3_orig))

  size_trim_expected <- fRagmentomics:::get_fragment_size(rs5_trim, rs3_trim)

  expect_identical(res_trim$Fragment_Size, as.integer(size_trim_expected))
  expect_lte(res_trim$Fragment_Size, res_no_trim$Fragment_Size)
})
