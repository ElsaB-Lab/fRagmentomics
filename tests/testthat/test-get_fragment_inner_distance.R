test_that("test get fragment inner distance", {
  # example 1
  r_pos1_1 <- 100
  r_cigar1_1 <- "2M6I2M"
  read_length1_1 <- 10
  r_pos2_1 <- 106
  r_cigar2_1 <- "10M"
  read_length2_1 <- 10
  inner_distance_1 <- get_fragment_inner_distance(r_pos1 = r_pos1_1, r_cigar1 = r_cigar1_1, read_length1 = read_length1_1, r_pos2 = r_pos2_1, r_cigar2 = r_cigar2_1, read_length2 = read_length2_1)

  expect_setequal(inner_distance_1, 2)

  # example 2
  r_pos1_2 <- 100
  r_cigar1_2 <- "1S2M6I2M1S"
  read_length1_2 <- 12
  r_pos2_2 <- 106
  r_cigar2_2 <- "12M"
  read_length2_2 <- 12
  inner_distance_2 <- get_fragment_inner_distance(r_pos1 = r_pos1_2, r_cigar1 = r_cigar1_2, read_length1 = read_length1_2, r_pos2 = r_pos2_2, r_cigar2 = r_cigar2_2, read_length2 = read_length2_2)

  expect_setequal(inner_distance_2, 1)


  # example 3
  r_pos1_3 <- 100
  r_cigar1_3 <- "1S2M6I3M"
  read_length1_3 <- 12
  r_pos2_3 <- 106
  r_cigar2_3 <- "1S11M"
  read_length2_3 <- 12
  inner_distance_3 <- get_fragment_inner_distance(r_pos1 = r_pos1_3, r_cigar1 = r_cigar1_3, read_length1 = read_length1_3, r_pos2 = r_pos2_3, r_cigar2 = r_cigar2_3, read_length2 = read_length2_3)

  expect_setequal(inner_distance_3, 0)


  # example 4
  r_pos1_4 <- 7578167
  r_cigar1_4 <- "125M1D14M5S"
  read_length1_4 <- 144
  r_pos2_4 <- 7578196
  r_cigar2_4 <- "17S96M1D31M"
  read_length2_4 <- 144
  inner_distance_4 <- get_fragment_inner_distance(r_pos1 = r_pos1_4, r_cigar1 = r_cigar1_4, read_length1 = read_length1_4, r_pos2 = r_pos2_4, r_cigar2 = r_cigar2_4, read_length2 = read_length2_4)

  expect_setequal(inner_distance_4, -133)


  # example 5
  r_pos1_5 <- 7578173
  r_cigar1_5 <- "40M1I103M"
  read_length1_5 <- 144
  r_pos2_5 <- 7578108
  r_cigar2_5 <- "105M1I38M"
  read_length2_5 <- 144
  inner_distance_5 <- get_fragment_inner_distance(r_pos1 = r_pos1_5, r_cigar1 = r_cigar1_5, read_length1 = read_length1_5, r_pos2 = r_pos2_5, r_cigar2 = r_cigar2_5, read_length2 = read_length2_5)

  expect_setequal(inner_distance_5, -78)


  # example 6
  r_pos1_6 <- 7578112
  r_cigar1_6 <- "98M46S"
  read_length1_6 <- 144
  r_pos2_6 <- 7578116
  r_cigar2_6 <- "33S94M17S"
  read_length2_6 <- 144
  inner_distance_6 <- get_fragment_inner_distance(r_pos1 = r_pos1_6, r_cigar1 = r_cigar1_6, read_length1 = read_length1_6, r_pos2 = r_pos2_6, r_cigar2 = r_cigar2_6, read_length2 = read_length2_6)

  expect_setequal(inner_distance_6, -173)
})
