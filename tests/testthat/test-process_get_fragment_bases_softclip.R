test_that("process_get_fragment_bases_softclip", {
  # Case 1: report_softclip == TRUE, both cigar strings have soft clips.
  res1 <- process_get_fragment_bases_softclip("45S100M", "100M30S")
  expect_equal(res1$nb_softclip_5p, 45)
  expect_equal(res1$nb_softclip_3p, 30)

  # Case 2: report_softclip == TRUE, cigar_5p has no soft clip while cigar_3p has one.
  res2 <- process_get_fragment_bases_softclip("100M", "20S100M30S")
  expect_equal(res2$nb_softclip_5p, 0)
  expect_equal(res2$nb_softclip_3p, 30)

  # Case 3: report_softclip == TRUE, cigar_5p has a soft clip while cigar_3p does not.
  res3 <- process_get_fragment_bases_softclip("45S100M", "100M")
  expect_equal(res3$nb_softclip_5p, 45)
  expect_equal(res3$nb_softclip_3p, 0)

  # Case 4: report_softclip == TRUE, neither cigar string contains a soft clip.
  res4 <- process_get_fragment_bases_softclip("100M", "100M")
  expect_equal(res4$nb_softclip_5p, 0)
  expect_equal(res4$nb_softclip_3p, 0)
})
