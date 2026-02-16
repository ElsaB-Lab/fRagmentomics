library(testthat)

# --- Helper to generate synthetic read dataframes ---
# FLAG defaults: 83L (first-in-pair, reverse strand) + 163L (second-in-pair, forward strand)
# This represents a properly oriented, valid R1/R2 pair.
make_reads <- function(n = 2, rname = "chr1", rnext = "=", pos = 100, tlen = 167,
                       flag = NULL) {
    if (is.null(flag)) {
        if (n == 0) {
            flag <- integer(0)
        } else if (n == 1) {
            flag <- 83L
        } else if (n == 2) {
            flag <- c(83L, 163L)
        } else {
            flag <- rep(83L, n)
        }
    }
    data.frame(
        RNAME = rep(rname, n),
        RNEXT = rep(rnext, n),
        POS = rep(pos, n),
        TLEN = rep(tlen, n),
        FLAG = as.integer(flag),
        stringsAsFactors = FALSE
    )
}

test_that("Case A: Returns correct diagnosis for exactly 1 read", {
    # 1. Mate is unmapped (RNEXT is *)
    df_unmapped <- make_reads(n = 1, rnext = "*")
    expect_match(
        process_fragment_reads_qc(df_unmapped, "chr1"),
        "Fragment has 1 read - Mate is unmapped"
    )

    # 1b. Mate is unmapped (RNEXT is NA)
    df_unmapped_na <- make_reads(n = 1, rnext = NA)
    expect_match(
        process_fragment_reads_qc(df_unmapped_na, "chr1"),
        "Fragment has 1 read - Mate is unmapped"
    )

    # 2. Translocation (Mate on different chromosome)
    df_transloc <- make_reads(n = 1, rnext = "chr2")
    expect_match(
        process_fragment_reads_qc(df_transloc, "chr1"),
        "Fragment has 1 read - Mate maps to a different chromosome: chr2"
    )

    # 3. Far Away (Mate on same chr but outside window) - with TLEN
    df_far <- make_reads(n = 1, rnext = "=", tlen = 5000)
    expect_match(
        process_fragment_reads_qc(df_far, "chr1"),
        "Fragment has 1 read - Mate maps outside loaded region. TLEN: 5000"
    )

    # 3b. Far Away - with NA TLEN (Covers the internal check for is.na(TLEN))
    df_far_na <- make_reads(n = 1, rnext = "=", tlen = NA)
    expect_match(
        process_fragment_reads_qc(df_far_na, "chr1"),
        "Fragment has 1 read - Mate maps outside loaded region. TLEN: NA"
    )
})

test_that("Case A: Returns correct diagnosis for abnormal counts != 1 and != 2", {
    # 0 Reads
    df_empty <- make_reads(n = 0)
    expect_equal(
        process_fragment_reads_qc(df_empty, "chr1"),
        "Fragment has 0 read(s)"
    )

    # 3 Reads
    df_triplet <- make_reads(n = 3)
    expect_equal(
        process_fragment_reads_qc(df_triplet, "chr1"),
        "Fragment has 3 read(s)"
    )
})

test_that("Case B: Returns empty string (OK) for valid 2-read fragments", {
    df_valid <- make_reads(n = 2, rname = "chr1", rnext = "=", pos = 100)
    expect_equal(
        process_fragment_reads_qc(df_valid, "chr1"),
        ""
    )
})

test_that("Case B: Returns correct diagnosis for 2 reads consistency checks", {
    # Test 1: Wrong Chromosome (Reads loaded but RNAME != expected chr)
    df_wrong_chr <- make_reads(n = 2, rname = "chr2")
    expect_match(
        process_fragment_reads_qc(df_wrong_chr, "chr1"),
        "Read\\(s\\) found on a chromosome other than chr1"
    )

    # Test 2: Missing Mate Info (RNEXT is *)
    df_no_mate_info <- make_reads(n = 2, rnext = "*")
    expect_match(
        process_fragment_reads_qc(df_no_mate_info, "chr1"),
        "Mate chromosome info is not available"
    )

    # Test 3: Internal Consistency (Reads present, but RNEXT points elsewhere)
    # Note: RNEXT is "chr2", so mate_chr != read$RNAME ("chr1")
    df_inconsistent <- make_reads(n = 2, rname = "chr1", rnext = "chr2")
    expect_match(
        process_fragment_reads_qc(df_inconsistent, "chr1"),
        "- Mate maps to a different chromosome: chr2"
    )

    # Test 4: Mapping Status (POS is 0)
    df_unmapped_pos <- make_reads(n = 2, pos = 0)
    expect_match(
        process_fragment_reads_qc(df_unmapped_pos, "chr1"),
        "One or both reads are unmapped \\(POS=0 or NA\\)"
    )

    # Test 4b: Mapping Status (POS is NA)
    df_na_pos <- make_reads(n = 2, pos = NA)
    expect_match(
        process_fragment_reads_qc(df_na_pos, "chr1"),
        "One or both reads are unmapped \\(POS=0 or NA\\)"
    )
})

test_that("Case B: Handles multiple errors simultaneously", {
    # Scenario: Wrong chromosome AND Unmapped POS
    df_multi <- data.frame(
        RNAME = c("chr2", "chr2"),
        RNEXT = c("=", "="),
        POS = c(0L, 0L),
        TLEN = c(100L, -100L),
        FLAG = c(83L, 163L),
        stringsAsFactors = FALSE
    )

    res <- process_fragment_reads_qc(df_multi, "chr1")

    expect_match(res, "found on a chromosome other than chr1")
    expect_match(res, "&") # Checks that errors are concatenated
    expect_match(res, "One or both reads are unmapped")
})
