library(testthat)

test_that("test identify_reads function:", {
  # Create test input: one read is first mate (flag 99), other is second mate (flag 147)
  # Flag 99 = isPaired + isProperPair + isFirstMateRead + isMateMinusStrand
  # Flag 147 = isPaired + isProperPair + isSecondMateRead + isMinusStrand
  df1 <- data.frame(
    qname = c("readA/1", "readA/2"),
    FLAG = c(99, 147),
    stringsAsFactors = FALSE
  )

  result1 <- identify_reads(df1)

  expect_equal(result1$read1$qname, "readA/1")
  expect_equal(result1$read2$qname, "readA/2")

  # identify_reads works with read1 and read2 in reverse order
  df2 <- data.frame(
    qname = c("readA/2", "readA/1"),
    FLAG = c(147, 99),
    stringsAsFactors = FALSE
  )

  result2 <- identify_reads(df2)

  expect_equal(result2$read1$qname, "readA/1")
  expect_equal(result2$read2$qname, "readA/2")

  # identify_reads throws error on invalid flags
  df3 <- data.frame(
    qname = c("readA/1", "readA/2"),
    FLAG = c(99, 99), # Both marked as first mate
    stringsAsFactors = FALSE
  )

  expect_error(identify_reads(df3), "Invalid read flags")
})
