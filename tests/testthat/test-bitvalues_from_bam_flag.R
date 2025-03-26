test_that("bitvalues from bam flag", {
  # example 1
  flag <- as.integer(163)
  bitnames <- c("isFirstMateRead", "isSecondMateRead")
  flag_1 <- bitvalues_from_bam_flag(flag, bitnames)

  expect_true(flag_1[, "isFirstMateRead"] == 0)
  expect_true(flag_1[, "isSecondMateRead"] == 1)

  # example 2
  flag <- as.integer(83)
  bitnames <- c("isFirstMateRead", "isSecondMateRead")
  flag_2 <- bitvalues_from_bam_flag(flag, bitnames)

  expect_true(flag_2[, "isFirstMateRead"] == 1)
  expect_true(flag_2[, "isSecondMateRead"] == 0)
})
