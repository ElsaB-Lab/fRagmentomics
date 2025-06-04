test_that("get_fragment_bases_5p_3p", {
  # Test 1: Standard case with n_bases less than sequence length
  res1 <- get_fragment_bases_5p_3p(
    n_bases = 4,
    seq_5p = "ACGTAAAAAAAAAAAAAAAAAAAAA",
    seq_3p = "AAAAAAAAAAAAAAAAAAAAATGCA",
    qual_5p = "IIII---------------------",
    qual_3p = "---------------------JJJJ"
  )
  expect_equal(res1$fragment_bases_5p, "ACGT")
  expect_equal(res1$fragment_bases_3p, "TGCA")
  expect_equal(res1$fragment_basqs_5p, "IIII")
  expect_equal(res1$fragment_basqs_3p, "JJJJ")

  # Test 2: n_bases equal bigger than the full length of the sequences
  expect_warning(
    {
      res2 <- get_fragment_bases_5p_3p(
        n_bases = 30,
        seq_5p = "ACGTAAAAAAAAAAAAAAAAAAAAA",
        seq_3p = "AAAAAAAAAAAAAAAAAAAAATGCA",
        qual_5p = "IIII---------------------",
        qual_3p = "---------------------JJJJ"
      )
    },
    regexp = "Number of bases in 5p and 3p bigger than the fragment length"
  )

  expect_equal(res2$fragment_bases_5p, "ACGTAAAAAAAAAAAAAAAAAAAAA")
  expect_equal(res2$fragment_bases_3p, "AAAAAAAAAAAAAAAAAAAAATGCA")
  expect_equal(res2$fragment_basqs_5p, "IIII---------------------")
  expect_equal(res2$fragment_basqs_3p, "---------------------JJJJ")

  # Test 3: n_bases = 0 should return empty strings for all fragments
  res3 <- get_fragment_bases_5p_3p(
    n_bases = 0,
    seq_5p = "ACGT",
    seq_3p = "TGCATGCA",
    qual_5p = "IIII",
    qual_3p = "JJJJ"
  )
  expect_equal(res3$fragment_bases_5p, "")
  expect_equal(res3$fragment_bases_3p, "")
  expect_equal(res3$fragment_basqs_5p, "")
  expect_equal(res3$fragment_basqs_3p, "")

  # Test 4: Standard case with n_bases == 1
  res1 <- get_fragment_bases_5p_3p(
    n_bases = 1,
    seq_5p = "CAAAAAAAAAAAAAAAAAAAAAAAA",
    seq_3p = "AAAAAAAAAAAAAAAAAAAAAAAAC",
    qual_5p = "I------------------------",
    qual_3p = "------------------------J"
  )
  expect_equal(res1$fragment_bases_5p, "C")
  expect_equal(res1$fragment_bases_3p, "C")
  expect_equal(res1$fragment_basqs_5p, "I")
  expect_equal(res1$fragment_basqs_3p, "J")
})
