# Helper: build a minimal SAM-like data.frame for process_fragment_reads_qc
.make_qc_df <- function(n = 2,
                        qnames = paste0("frag", seq_len(n)),
                        flags = NULL,
                        rnames = rep("chr1", n),
                        pos = seq(100L, by = 100L, length.out = n),
                        rnext = rep("=", n),
                        pnext = rev(pos),
                        tlen = rep(200L, n)) {
  # Default flags: proper pair (83 + 163)
  if (is.null(flags)) {
    if (n == 1) {
      flags <- 83L
    } else if (n == 2) {
      flags <- c(83L, 163L)
    } else {
      flags <- rep(83L, n)
    }
  }
  data.frame(
    QNAME = rep("fragA", n),
    FLAG = as.integer(flags),
    RNAME = rnames,
    POS = as.integer(pos),
    TLEN = as.integer(tlen),
    MAPQ = rep(60L, n),
    CIGAR = rep("10M", n),
    RNEXT = rnext,
    PNEXT = as.integer(pnext),
    SEQ = rep("ACGTACGTAA", n),
    QUAL = rep("FFFFFFFFFF", n),
    stringsAsFactors = FALSE
  )
}


# ---- Case: 2 reads, all OK --------------------------------------------------

test_that("2 well-formed reads return empty QC string", {
  # FLAG 83 (first, reverse) + FLAG 163 (second, forward) = proper pair
  df <- .make_qc_df(n = 2, flags = c(83L, 163L))
  result <- fRagmentomics:::process_fragment_reads_qc(df, "chr1")
  expect_equal(result, "")
})


# ---- Case A: Abnormal read count --------------------------------------------

test_that("0 reads returns appropriate message", {
  df <- .make_qc_df(n = 2)[FALSE, , drop = FALSE] # empty df with right columns
  result <- fRagmentomics:::process_fragment_reads_qc(df, "chr1")
  expect_match(result, "Fragment has 0 read")
})

test_that("1 read with unmapped mate returns diagnosis", {
  # FLAG 83: isMinusStrand=1, isMateMinusStrand=0 (well oriented)
  df <- .make_qc_df(n = 1, flags = 83L, rnext = "*")
  result <- fRagmentomics:::process_fragment_reads_qc(df, "chr1")
  expect_match(result, "Fragment has 1 read")
  expect_match(result, "Mate is unmapped")
})

test_that("1 read with mate on different chromosome returns translocation", {
  # FLAG 83: well oriented
  df <- .make_qc_df(n = 1, flags = 83L, rnext = "chr2")
  result <- fRagmentomics:::process_fragment_reads_qc(df, "chr1")
  expect_match(result, "Fragment has 1 read")
  expect_match(result, "Mate maps to a different chromosome: chr2")
})

test_that("1 read with mate far away returns TLEN diagnosis", {
  # FLAG 83: well oriented, RNEXT = '=' (same chr), mate just outside window
  df <- .make_qc_df(n = 1, flags = 83L, rnext = "=", tlen = 5000L)
  result <- fRagmentomics:::process_fragment_reads_qc(df, "chr1")
  expect_match(result, "Fragment has 1 read")
  expect_match(result, "Mate maps outside loaded region")
  expect_match(result, "5000")
})

test_that("3 reads returns appropriate message", {
  df <- .make_qc_df(n = 3, flags = c(83L, 163L, 83L))
  result <- fRagmentomics:::process_fragment_reads_qc(df, "chr1")
  expect_match(result, "Fragment has 3 read")
})


# ---- Case B: 2-read consistency checks --------------------------------------

test_that("2 reads on wrong chromosome are flagged", {
  df <- .make_qc_df(n = 2, flags = c(83L, 163L), rnames = c("chr2", "chr2"))
  result <- fRagmentomics:::process_fragment_reads_qc(df, "chr1")
  expect_match(result, "Read\\(s\\) found on a chromosome other than chr1")
})

test_that("2 reads with missing mate info (RNEXT = *) are flagged", {
  df <- .make_qc_df(n = 2, flags = c(83L, 163L), rnext = c("*", "="))
  result <- fRagmentomics:::process_fragment_reads_qc(df, "chr1")
  expect_match(result, "Mate chromosome info is not available")
})

test_that("2 reads with inconsistent RNEXT are flagged", {
  df <- .make_qc_df(n = 2, flags = c(83L, 163L), rnext = c("chr2", "="))
  result <- fRagmentomics:::process_fragment_reads_qc(df, "chr1")
  expect_match(result, "Mate maps to a different chromosome: chr2")
})

test_that("2 reads with unmapped position (POS=0) are flagged", {
  df <- .make_qc_df(n = 2, flags = c(83L, 163L), pos = c(0L, 200L))
  result <- fRagmentomics:::process_fragment_reads_qc(df, "chr1")
  expect_match(result, "One or both reads are unmapped")
})

test_that("2 reads with invalid R1/R2 pairing (both first mate) are flagged", {
  # FLAG 83 = first in pair, reverse.  Use two first-mate reads.
  df <- .make_qc_df(n = 2, flags = c(83L, 83L))
  result <- fRagmentomics:::process_fragment_reads_qc(df, "chr1")
  expect_match(result, "Fragment is not a valid R1/R2 pair")
})

test_that("2 reads badly oriented (both same strand) are flagged", {
  # FLAG 67 = paired(1) + proper(2) + first(64) = forward, mate forward
  # FLAG 131 = paired(1) + proper(2) + second(128) = forward, mate forward
  # Both isMinusStrand=0 -> sum=0 != 1 -> badly oriented
  df <- .make_qc_df(n = 2, flags = c(67L, 131L))
  result <- fRagmentomics:::process_fragment_reads_qc(df, "chr1")
  expect_match(result, "Badly oriented reads")
})


# ---- Case C: 1-read orientation check ---------------------------------------

test_that("1 read with badly oriented mate is flagged", {
  # FLAG 67 = paired(1) + proper(2) + first(64) = forward, mate forward
  # isMinusStrand=0, isMateMinusStrand=0 -> equal -> badly oriented
  df <- .make_qc_df(n = 1, flags = 67L, rnext = "=", tlen = 5000L)
  result <- fRagmentomics:::process_fragment_reads_qc(df, "chr1")
  expect_match(result, "Badly oriented reads")
})

test_that("1 read with well oriented mate does not flag orientation", {
  # FLAG 83 = first in pair, reverse, mate forward
  # isMinusStrand=1, isMateMinusStrand=0 -> different -> OK
  df <- .make_qc_df(n = 1, flags = 83L, rnext = "=", tlen = 5000L)
  result <- fRagmentomics:::process_fragment_reads_qc(df, "chr1")
  expect_no_match(result, "Badly oriented reads")
})


# ---- Additive messages -------------------------------------------------------

test_that("Multiple QC failures are combined with ' & '", {
  # 2 reads: wrong chr + both same strand + invalid R1/R2
  # FLAG 67 = paired + proper + first + forward + mate forward
  # -> isFirstMateRead=1 for both -> sum=2 -> invalid pair
  # -> isMinusStrand=0 for both -> sum=0 -> badly oriented
  df <- .make_qc_df(
    n = 2,
    flags = c(67L, 67L),
    rnames = c("chr2", "chr2")
  )
  result <- fRagmentomics:::process_fragment_reads_qc(df, "chr1")

  # All three issues should be present, separated by " & "
  expect_match(result, "Read\\(s\\) found on a chromosome other than chr1")
  expect_match(result, "Fragment is not a valid R1/R2 pair")
  expect_match(result, "Badly oriented reads")
  expect_match(result, " & ")
})

test_that("1 read: far away + badly oriented are combined", {
  # FLAG 67: forward, mate forward -> isMinusStrand==isMateMinusStrand -> badly oriented
  # RNEXT = '=' -> same chr -> mate far away (1 read case)
  df <- .make_qc_df(n = 1, flags = 67L, rnext = "=", tlen = 8000L)
  result <- fRagmentomics:::process_fragment_reads_qc(df, "chr1")

  expect_match(result, "Fragment has 1 read")
  expect_match(result, "Mate maps outside loaded region")
  expect_match(result, "Badly oriented reads")
  expect_match(result, " & ")
})
