test_that("process_get_fragment_bases", {
  # Test 1: Standard case with nbr_bases less than sequence length
  res1 <- process_get_fragment_bases(
    nbr_bases = 4,
    query_5p = "ACGTAAAAAAAAAAAAAAAAAAAAA",
    query_3p = "AAAAAAAAAAAAAAAAAAAAATGCA",
    qual_5p = "IIII---------------------",
    qual_3p = "---------------------JJJJ"
  )
  expect_equal(res1$fragment_bases_5p, "ACGT")
  expect_equal(res1$fragment_bases_3p, "TGCA")
  expect_equal(res1$fragment_Qbases_5p, "IIII")
  expect_equal(res1$fragment_Qbases_3p, "JJJJ")

  # Test 2: nbr_bases equal bigger than the full length of the sequences
  expect_warning(
    {
      res2 <- process_get_fragment_bases(
        nbr_bases = 30,
        query_5p = "ACGTAAAAAAAAAAAAAAAAAAAAA",
        query_3p = "AAAAAAAAAAAAAAAAAAAAATGCA",
        qual_5p = "IIII---------------------",
        qual_3p = "---------------------JJJJ"
      )
    },
    regexp = "Number of bases in 5p and 3p bigger than the fragment length"
  )

  expect_equal(res2$fragment_bases_5p, "ACGTAAAAAAAAAAAAAAAAAAAAA")
  expect_equal(res2$fragment_bases_3p, "AAAAAAAAAAAAAAAAAAAAATGCA")
  expect_equal(res2$fragment_Qbases_5p, "IIII---------------------")
  expect_equal(res2$fragment_Qbases_3p, "---------------------JJJJ")

  # Test 3: nbr_bases = 0 should return empty strings for all fragments
  res3 <- process_get_fragment_bases(
    nbr_bases = 0,
    query_5p = "ACGT",
    query_3p = "TGCATGCA",
    qual_5p = "IIII",
    qual_3p = "JJJJ"
  )
  expect_equal(res3$fragment_bases_5p, "")
  expect_equal(res3$fragment_bases_3p, "")
  expect_equal(res3$fragment_Qbases_5p, "")
  expect_equal(res3$fragment_Qbases_3p, "")

  # Test 4: Standard case with nbr_bases == 1
  res1 <- process_get_fragment_bases(
    nbr_bases = 1,
    query_5p = "CAAAAAAAAAAAAAAAAAAAAAAAA",
    query_3p = "AAAAAAAAAAAAAAAAAAAAAAAAC",
    qual_5p = "I------------------------",
    qual_3p = "------------------------J"
  )
  expect_equal(res1$fragment_bases_5p, "C")
  expect_equal(res1$fragment_bases_3p, "C")
  expect_equal(res1$fragment_Qbases_5p, "I")
  expect_equal(res1$fragment_Qbases_3p, "J")
})
