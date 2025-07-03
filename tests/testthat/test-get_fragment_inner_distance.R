test_that("test get fragment inner distance", {
  # example 1
  r_pos_5p_1 <- 100
  r_cigar_5p_1 <- "2M6I2M"
  r_pos_3p_1 <- 106
  r_cigar_3p_1 <- "10M"
  inner_distance_1 <- get_fragment_inner_distance(r_pos_5p = r_pos_5p_1, r_cigar_5p = r_cigar_5p_1, r_pos_3p = r_pos_3p_1, r_cigar_3p = r_cigar_3p_1)

  expect_setequal(inner_distance_1, 2)

  # example 2
  r_pos_5p_2 <- 100
  r_cigar_5p_2 <- "1S2M6I2M1S"
  r_pos_3p_2 <- 106
  r_cigar_3p_2 <- "12M"
  inner_distance_2 <- get_fragment_inner_distance(r_pos_5p = r_pos_5p_2, r_cigar_5p = r_cigar_5p_2, r_pos_3p = r_pos_3p_2, r_cigar_3p = r_cigar_3p_2)

  expect_setequal(inner_distance_2, 1)


  # example 3
  r_pos_5p_3 <- 100
  r_cigar_5p_3 <- "1S2M6I3M"
  r_pos_3p_3 <- 106
  r_cigar_3p_3 <- "1S11M"
  inner_distance_3 <- get_fragment_inner_distance(r_pos_5p = r_pos_5p_3, r_cigar_5p = r_cigar_5p_3, r_pos_3p = r_pos_3p_3, r_cigar_3p = r_cigar_3p_3)

  expect_setequal(inner_distance_3, 0)


  # example 4
  r_pos_5p_4 <- 7578167
  r_cigar_5p_4 <- "125M1D14M5S"
  r_pos_3p_4 <- 7578196
  r_cigar_3p_4 <- "17S96M1D31M"
  inner_distance_4 <- get_fragment_inner_distance(r_pos_5p = r_pos_5p_4, r_cigar_5p = r_cigar_5p_4, r_pos_3p = r_pos_3p_4, r_cigar_3p = r_cigar_3p_4)

  expect_setequal(inner_distance_4, -133)


  # example 5
  r_pos_3p_5 <- 7578173
  r_cigar_3p_5 <- "40M1I103M"
  r_pos_5p_5 <- 7578108
  r_cigar_5p_5 <- "105M1I38M"
  inner_distance_5 <- get_fragment_inner_distance(r_pos_5p = r_pos_5p_5, r_cigar_5p = r_cigar_5p_5, r_pos_3p = r_pos_3p_5, r_cigar_3p = r_cigar_3p_5)

  expect_setequal(inner_distance_5, -78)


  # example 6
  r_pos_5p_6 <- 7578112
  r_cigar_5p_6 <- "98M46S"
  r_pos_3p_6 <- 7578116
  r_cigar_3p_6 <- "33S94M17S"
  inner_distance_6 <- get_fragment_inner_distance(r_pos_5p = r_pos_5p_6, r_cigar_5p = r_cigar_5p_6, r_pos_3p = r_pos_3p_6, r_cigar_3p = r_cigar_3p_6)

  expect_setequal(inner_distance_6, -173)
})
