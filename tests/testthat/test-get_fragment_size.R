test_that("test get fragment size", {
  # example 1
  r_cigar_5p_1 <- "2M6I2M"
  r_cigar_3p_1 <- "10M"
  read_stats_5p_1 <- list(
    POS = 100, CIGAR = r_cigar_5p_1, read_length = calculate_read_length(r_cigar_5p_1)
  )
  read_stats_3p_1 <- list(
    POS = 106, CIGAR = r_cigar_3p_1, read_length = calculate_read_length(r_cigar_3p_1)
  )

  expect_equal(get_fragment_size(read_stats_5p_1, read_stats_3p_1), 22)

  # example 2
  r_cigar_5p_2 <- "1S2M6I2M1S"
  r_cigar_3p_2 <- "12M"

  read_stats_5p_2 <- list(
    POS = 100, CIGAR = r_cigar_5p_2, read_length = calculate_read_length(r_cigar_5p_2)
  )
  read_stats_3p_2 <- list(
    POS = 106, CIGAR = r_cigar_3p_2, read_length = calculate_read_length(r_cigar_3p_2)
  )

  expect_equal(get_fragment_size(read_stats_5p_2, read_stats_3p_2), 25)


  # example 3
  r_cigar_5p_3 <- "1S2M6I3M"
  r_cigar_3p_3 <- "1S11M"

  read_stats_5p_3 <- list(
    POS = 100, CIGAR = r_cigar_5p_3, read_length = calculate_read_length(r_cigar_5p_3)
  )
  read_stats_3p_3 <- list(
    POS = 106, CIGAR = r_cigar_3p_3, read_length = calculate_read_length(r_cigar_3p_3)
  )

  expect_equal(get_fragment_size(read_stats_5p_3, read_stats_3p_3), 24)


  # example 4
  r_cigar_5p_4 <- "125M1D14M5S"
  r_cigar_3p_4 <- "17S96M1D31M"

  read_stats_5p_4 <- list(
    POS = 7578167, CIGAR = r_cigar_5p_4, read_length = calculate_read_length(r_cigar_5p_4)
  )
  read_stats_3p_4 <- list(
    POS = 7578196, CIGAR = r_cigar_3p_4, read_length = calculate_read_length(r_cigar_3p_4)
  )

  expect_equal(get_fragment_size(read_stats_5p_4, read_stats_3p_4), 156)


  # example 5
  r_cigar_5p_5 <- "105M1I38M"
  r_cigar_3p_5 <- "40M1I103M"

  read_stats_5p_5 <- list(
    POS = 7578108, CIGAR = r_cigar_5p_5, read_length = calculate_read_length(r_cigar_5p_5)
  )
  read_stats_3p_5 <- list(
    POS = 7578173, CIGAR = r_cigar_3p_5, read_length = calculate_read_length(r_cigar_3p_5)
  )

  expect_equal(get_fragment_size(read_stats_5p_5, read_stats_3p_5), 209)


  # example 6
  r_cigar_5p_6 <- "98M46S"
  r_cigar_3p_6 <- "33S94M17S"

  read_stats_5p_6 <- list(
    POS = 7578112, CIGAR = r_cigar_5p_6, read_length = calculate_read_length(r_cigar_5p_6)
  )
  read_stats_3p_6 <- list(
    POS = 7578116, CIGAR = r_cigar_3p_6, read_length = calculate_read_length(r_cigar_3p_6)
  )

  expect_equal(get_fragment_size(read_stats_5p_6, read_stats_3p_6), 115)

  # Deletion
  # Del 1 - Normal case, deletion in overlap
  r_cigar_5p_d1 <- "4M4D4M"
  r_cigar_3p_d1 <- "1M4D7M"
  read_stats_5p_d1 <- list(
    POS = 1, CIGAR = r_cigar_5p_d1, read_length = calculate_read_length(r_cigar_5p_d1)
  )
  read_stats_3p_d1 <- list(
    POS = 4, CIGAR = r_cigar_3p_d1, read_length = calculate_read_length(r_cigar_3p_d1)
  )
  expect_equal(get_fragment_size(read_stats_5p_d1, read_stats_3p_d1), 11)

  # Del 2 - Normal case, deletion on 5p and 3p no overlap
  r_cigar_5p_d2 <- "4M3D1M"
  r_cigar_3p_d2 <- "1M3D4M"
  read_stats_5p_d2 <- list(
    POS = 1, CIGAR = r_cigar_5p_d2, read_length = calculate_read_length(r_cigar_5p_d2)
  )
  read_stats_3p_d2 <- list(
    POS = 12, CIGAR = r_cigar_3p_d2, read_length = calculate_read_length(r_cigar_3p_d2)
  )
  expect_equal(get_fragment_size(read_stats_5p_d2, read_stats_3p_d2), 13)

  # Del 3 - Error sequencing 5p
  r_cigar_5p_d3 <- "7M6D1M"
  r_cigar_3p_d3 <- "8M"
  read_stats_5p_d3 <- list(
    POS = 1, CIGAR = r_cigar_5p_d3, read_length = calculate_read_length(r_cigar_5p_d3)
  )
  read_stats_3p_d3 <- list(
    POS = 7, CIGAR = r_cigar_3p_d3, read_length = calculate_read_length(r_cigar_3p_d3)
  )
  expect_equal(get_fragment_size(read_stats_5p_d3, read_stats_3p_d3), 8)

  # Del 4 - Error sequencing overlap 5p
  r_cigar_5p_d4 <- "7M6D1M"
  r_cigar_3p_d4 <- "8M"
  read_stats_5p_d4 <- list(
    POS = 1, CIGAR = r_cigar_5p_d4, read_length = calculate_read_length(r_cigar_5p_d4)
  )
  read_stats_3p_d4 <- list(
    POS = 11, CIGAR = r_cigar_3p_d4, read_length = calculate_read_length(r_cigar_3p_d4)
  )
  expect_equal(get_fragment_size(read_stats_5p_d4, read_stats_3p_d4), 12)

  # Del 5 - Error sequencing 3p
  r_cigar_5p_d5 <- "8M"
  r_cigar_3p_d5 <- "1M5D7M"
  read_stats_5p_d5 <- list(
    POS = 1, CIGAR = r_cigar_5p_d5, read_length = calculate_read_length(r_cigar_5p_d5)
  )
  read_stats_3p_d5 <- list(
    POS = 3, CIGAR = r_cigar_3p_d5, read_length = calculate_read_length(r_cigar_3p_d5)
  )
  expect_equal(get_fragment_size(read_stats_5p_d5, read_stats_3p_d5), 10)

  # Del 6 - Error sequencing overlap 3p
  r_cigar_5p_d6 <- "8M"
  r_cigar_3p_d6 <- "1M5D7M"
  read_stats_5p_d6 <- list(
    POS = 1, CIGAR = r_cigar_5p_d6, read_length = calculate_read_length(r_cigar_5p_d6)
  )
  read_stats_3p_d6 <- list(
    POS = 6, CIGAR = r_cigar_3p_d6, read_length = calculate_read_length(r_cigar_3p_d6)
  )
  expect_equal(get_fragment_size(read_stats_5p_d6, read_stats_3p_d6), 13)

  # Del 7 - Error sequencing overlap and soft clipping
  r_cigar_5p_d7 <- "1M1D5M"
  r_cigar_3p_d7 <- "1S1M1D6M"
  read_stats_5p_d7 <- list(
    POS = 1, CIGAR = r_cigar_5p_d7, read_length = calculate_read_length(r_cigar_5p_d7)
  )
  read_stats_3p_d7 <- list(
    POS = 3, CIGAR = r_cigar_3p_d7, read_length = calculate_read_length(r_cigar_3p_d7)
  )
  expect_equal(get_fragment_size(read_stats_5p_d7, read_stats_3p_d7), 8)

  # Del 8 - Discordant case
  r_cigar_5p_d7 <- "4M4D4M"
  r_cigar_3p_d7 <- "1M6D7M"
  read_stats_5p_d7 <- list(
    POS = 1, CIGAR = r_cigar_5p_d7, read_length = calculate_read_length(r_cigar_5p_d7)
  )
  read_stats_3p_d7 <- list(
    POS = 4, CIGAR = r_cigar_3p_d7, read_length = calculate_read_length(r_cigar_3p_d7)
  )
  expect_equal(get_fragment_size(read_stats_5p_d7, read_stats_3p_d7), 11)

  # Del 9 - Several deletion in the overlapping section
  r_cigar_5p_d7 <- "1M1D1M1D4M"
  r_cigar_3p_d7 <- "1M1D1M1D4M"
  read_stats_5p_d7 <- list(
    POS = 1, CIGAR = r_cigar_5p_d7, read_length = calculate_read_length(r_cigar_5p_d7)
  )
  read_stats_3p_d7 <- list(
    POS = 1, CIGAR = r_cigar_3p_d7, read_length = calculate_read_length(r_cigar_3p_d7)
  )
  expect_equal(get_fragment_size(read_stats_5p_d7, read_stats_3p_d7), 6)

  # Del 10 - Several deletion in the overlapping section
  r_cigar_5p_d7 <- "1M1D1M1D4M"
  r_cigar_3p_d7 <- "1M1D1M1D4M"
  read_stats_5p_d7 <- list(
    POS = 1, CIGAR = r_cigar_5p_d7, read_length = calculate_read_length(r_cigar_5p_d7)
  )
  read_stats_3p_d7 <- list(
    POS = 2, CIGAR = r_cigar_3p_d7, read_length = calculate_read_length(r_cigar_3p_d7)
  )
  expect_equal(get_fragment_size(read_stats_5p_d7, read_stats_3p_d7), 5)

  # Insertion
  # Ins 1 - Normal case, insertion in overlap
  r_cigar_5p_i1 <- "4M4I4M"
  r_cigar_3p_i1 <- "1M4I7M"
  read_stats_5p_i1 <- list(
    POS = 1, CIGAR = r_cigar_5p_i1, read_length = calculate_read_length(r_cigar_5p_i1)
  )
  read_stats_3p_i1 <- list(
    POS = 4, CIGAR = r_cigar_3p_i1, read_length = calculate_read_length(r_cigar_3p_i1)
  )
  expect_equal(get_fragment_size(read_stats_5p_i1, read_stats_3p_i1), 15)

  # Ins 2 - Normal case, insertion on 5p and 3p no overlap
  r_cigar_5p_i2 <- "4M4I4M"
  r_cigar_3p_i2 <- "1M4I7M"
  read_stats_5p_i2 <- list(
    POS = 1, CIGAR = r_cigar_5p_i2, read_length = calculate_read_length(r_cigar_5p_i2)
  )
  read_stats_3p_i2 <- list(
    POS = 14, CIGAR = r_cigar_3p_i2, read_length = calculate_read_length(r_cigar_3p_i2)
  )
  expect_equal(get_fragment_size(read_stats_5p_i2, read_stats_3p_i2), 29)

  # Ins 3 - Overlap normal cases 3p - Are they existing ???
  r_cigar_5p_i3 <- "7M6I1M"
  r_cigar_3p_i3 <- "3I11M"
  read_stats_5p_i3 <- list(
    POS = 1, CIGAR = r_cigar_5p_i3, read_length = calculate_read_length(r_cigar_5p_i3)
  )
  read_stats_3p_i3 <- list(
    POS = 8, CIGAR = r_cigar_3p_i3, read_length = calculate_read_length(r_cigar_3p_i3)
  )
  expect_equal(get_fragment_size(read_stats_5p_i3, read_stats_3p_i3), 24)

  # Ins 4
  r_cigar_5p_i4 <- "10M4I"
  r_cigar_3p_i4 <- "2M6I6M"
  read_stats_5p_i4 <- list(
    POS = 1, CIGAR = r_cigar_5p_i4, read_length = calculate_read_length(r_cigar_5p_i4)
  )
  read_stats_3p_i4 <- list(
    POS = 9, CIGAR = r_cigar_3p_i4, read_length = calculate_read_length(r_cigar_3p_i4)
  )
  expect_equal(get_fragment_size(read_stats_5p_i4, read_stats_3p_i4), 22)

  # Ins 4 bis
  r_cigar_5p_i4 <- "10M4S"
  r_cigar_3p_i4 <- "2M6I6M"
  read_stats_5p_i4 <- list(
    POS = 1, CIGAR = r_cigar_5p_i4, read_length = calculate_read_length(r_cigar_5p_i4)
  )
  read_stats_3p_i4 <- list(
    POS = 9, CIGAR = r_cigar_3p_i4, read_length = calculate_read_length(r_cigar_3p_i4)
  )
  expect_equal(get_fragment_size(read_stats_5p_i4, read_stats_3p_i4), 22)

  # Ins 4 ter
  r_cigar_5p_i4 <- "10M4S"
  r_cigar_3p_i4 <- "2M3I9M"
  read_stats_5p_i4 <- list(
    POS = 1, CIGAR = r_cigar_5p_i4, read_length = calculate_read_length(r_cigar_5p_i4)
  )
  read_stats_3p_i4 <- list(
    POS = 9, CIGAR = r_cigar_3p_i4, read_length = calculate_read_length(r_cigar_3p_i4)
  )
  expect_equal(get_fragment_size(read_stats_5p_i4, read_stats_3p_i4), 22)

  # Ins 5 - Error sequencing 5p - Are they existing ???
  r_cigar_5p_i5 <- "4M4I4M"
  r_cigar_3p_i5 <- "12M"
  read_stats_5p_i5 <- list(
    POS = 1, CIGAR = r_cigar_5p_i5, read_length = calculate_read_length(r_cigar_5p_i5)
  )
  read_stats_3p_i5 <- list(
    POS = 4, CIGAR = r_cigar_3p_i5, read_length = calculate_read_length(r_cigar_3p_i5)
  )
  expect_equal(get_fragment_size(read_stats_5p_i5, read_stats_3p_i5), 19)

  # Ins 6 - Error sequencing 3p
  r_cigar_5p_i6 <- "12M"
  r_cigar_3p_i6 <- "1M4I7M"
  read_stats_5p_i6 <- list(
    POS = 1, CIGAR = r_cigar_5p_i6, read_length = calculate_read_length(r_cigar_5p_i6)
  )
  read_stats_3p_i6 <- list(
    POS = 4, CIGAR = r_cigar_3p_i6, read_length = calculate_read_length(r_cigar_3p_i6)
  )
  expect_equal(get_fragment_size(read_stats_5p_i6, read_stats_3p_i6), 15)

  # Ins 7 - Error sequencing overlap 5p partially
  r_cigar_5p_i7 <- "10M2I3S"
  r_cigar_3p_i7 <- "1M4I10M"
  read_stats_5p_i7 <- list(
    POS = 1, CIGAR = r_cigar_5p_i7, read_length = calculate_read_length(r_cigar_5p_i7)
  )
  read_stats_3p_i7 <- list(
    POS = 10, CIGAR = r_cigar_3p_i7, read_length = calculate_read_length(r_cigar_3p_i7)
  )
  expect_equal(get_fragment_size(read_stats_5p_i7, read_stats_3p_i7), 24)

  # Ins 8 - Complexe case with soft clipping, insertion and deletion
  r_cigar_5p_c8 <- "1M2I1M1D1M1S"
  r_cigar_3p_c8 <- "1S1M1D1M1D2M"
  read_stats_5p_c8 <- list(
    POS = 1, CIGAR = r_cigar_5p_c8, read_length = calculate_read_length(r_cigar_5p_c8)
  )
  read_stats_3p_c8 <- list(
    POS = 2, CIGAR = r_cigar_3p_c8, read_length = calculate_read_length(r_cigar_3p_c8)
  )
  expect_equal(get_fragment_size(read_stats_5p_c8, read_stats_3p_c8), 7)

  # Ins 9 - Error sequencing overlap 5p partially
  r_cigar_5p_i9 <- "7M6I1M"
  r_cigar_3p_i9 <- "3S11M"
  read_stats_5p_i9 <- list(
    POS = 1, CIGAR = r_cigar_5p_i9, read_length = calculate_read_length(r_cigar_5p_i9)
  )
  read_stats_3p_i9 <- list(
    POS = 8, CIGAR = r_cigar_3p_i9, read_length = calculate_read_length(r_cigar_3p_i9)
  )
  expect_equal(get_fragment_size(read_stats_5p_i9, read_stats_3p_i9), 24)

  # Ins 10 - Real case
  r_cigar_5p_i10 <- "138M6S"
  r_cigar_3p_i10 <- "37M35I72M"
  read_stats_5p_i10 <- list(
    POS = 868, CIGAR = r_cigar_5p_i10, read_length = calculate_read_length(r_cigar_5p_i10)
  )
  read_stats_3p_i10 <- list(
    POS = 936, CIGAR = r_cigar_3p_i10, read_length = calculate_read_length(r_cigar_3p_i10)
  )
  expect_equal(get_fragment_size(read_stats_5p_i10, read_stats_3p_i10), 212)

  # Soft Clipping
  # Soft clipping 1
  r_cigar_5p_s1 <- "3S4M"
  r_cigar_3p_s1 <- "4M3S"
  read_stats_5p_s1 <- list(
    POS = 4, CIGAR = r_cigar_5p_s1, read_length = calculate_read_length(r_cigar_5p_s1)
  )
  read_stats_3p_s1 <- list(
    POS = 5, CIGAR = r_cigar_3p_s1, read_length = calculate_read_length(r_cigar_3p_s1)
  )
  expect_equal(get_fragment_size(read_stats_5p_s1, read_stats_3p_s1), 11)

  # Soft clipping 2
  r_cigar_5p_s2 <- "3S4M"
  r_cigar_3p_s2 <- "4M3S"
  read_stats_5p_s2 <- list(
    POS = 4, CIGAR = r_cigar_5p_s2, read_length = calculate_read_length(r_cigar_5p_s2)
  )
  read_stats_3p_s2 <- list(
    POS = 11, CIGAR = r_cigar_3p_s2, read_length = calculate_read_length(r_cigar_3p_s2)
  )
  expect_equal(get_fragment_size(read_stats_5p_s2, read_stats_3p_s2), 17)

  # Soft clipping 3
  r_cigar_5p_s3 <- "2M2S"
  r_cigar_3p_s3 <- "2S2M"
  read_stats_5p_s3 <- list(
    POS = 5, CIGAR = r_cigar_5p_s3, read_length = calculate_read_length(r_cigar_5p_s3)
  )
  read_stats_3p_s3 <- list(
    POS = 5, CIGAR = r_cigar_3p_s3, read_length = calculate_read_length(r_cigar_3p_s3)
  )
  expect_equal(get_fragment_size(read_stats_5p_s3, read_stats_3p_s3), 2)

  # Soft clipping 4
  r_cigar_5p_s4 <- "2S4M"
  r_cigar_3p_s4 <- "4M2S"
  read_stats_5p_s4 <- list(
    POS = 5, CIGAR = r_cigar_5p_s4, read_length = calculate_read_length(r_cigar_5p_s4)
  )
  read_stats_3p_s4 <- list(
    POS = 3, CIGAR = r_cigar_3p_s4, read_length = calculate_read_length(r_cigar_3p_s4)
  )
  expect_equal(get_fragment_size(read_stats_5p_s4, read_stats_3p_s4), 6)

  # Soft clipping 5
  r_cigar_5p_s5 <- "6M"
  r_cigar_3p_s5 <- "4M2S"
  read_stats_5p_s5 <- list(
    POS = 3, CIGAR = r_cigar_5p_s5, read_length = calculate_read_length(r_cigar_5p_s5)
  )
  read_stats_3p_s5 <- list(
    POS = 3, CIGAR = r_cigar_3p_s5, read_length = calculate_read_length(r_cigar_3p_s5)
  )
  expect_equal(get_fragment_size(read_stats_5p_s5, read_stats_3p_s5), 6)

  # Soft clipping 6
  r_cigar_5p_s6 <- "6M"
  r_cigar_3p_s6 <- "2S4M"
  read_stats_5p_s6 <- list(
    POS = 3, CIGAR = r_cigar_5p_s6, read_length = calculate_read_length(r_cigar_5p_s6)
  )
  read_stats_3p_s6 <- list(
    POS = 5, CIGAR = r_cigar_3p_s6, read_length = calculate_read_length(r_cigar_3p_s6)
  )
  expect_equal(get_fragment_size(read_stats_5p_s6, read_stats_3p_s6), 6)

  # Soft clipping 7
  r_cigar_5p_s7 <- "6M"
  r_cigar_3p_s7 <- "2S4M"
  read_stats_5p_s7 <- list(
    POS = 3, CIGAR = r_cigar_5p_s7, read_length = calculate_read_length(r_cigar_5p_s7)
  )
  read_stats_3p_s7 <- list(
    POS = 2, CIGAR = r_cigar_3p_s7, read_length = calculate_read_length(r_cigar_3p_s7)
  )
  expect_equal(get_fragment_size(read_stats_5p_s7, read_stats_3p_s7), 3)

  # Soft clipping 8
  r_cigar_5p_s8 <- "6M"
  r_cigar_3p_s8 <- "4M2S"
  read_stats_5p_s8 <- list(
    POS = 5, CIGAR = r_cigar_5p_s8, read_length = calculate_read_length(r_cigar_5p_s8)
  )
  read_stats_3p_s8 <- list(
    POS = 3, CIGAR = r_cigar_3p_s8, read_length = calculate_read_length(r_cigar_3p_s8)
  )
  expect_equal(get_fragment_size(read_stats_5p_s8, read_stats_3p_s8), 4)

  # Soft clipping 9
  r_cigar_5p_s9 <- "2M3S"
  r_cigar_3p_s9 <- "2S3M"
  read_stats_5p_s9 <- list(
    POS = 5, CIGAR = r_cigar_5p_s9, read_length = calculate_read_length(r_cigar_5p_s9)
  )
  read_stats_3p_s9 <- list(
    POS = 5, CIGAR = r_cigar_3p_s9, read_length = calculate_read_length(r_cigar_3p_s9)
  )
  expect_equal(get_fragment_size(read_stats_5p_s9, read_stats_3p_s9), 3)

  # Soft clipping 10
  r_cigar_5p_s10 <- "3S4M2S"
  r_cigar_3p_s10 <- "3S4M3S"
  read_stats_5p_s10 <- list(
    POS = 4, CIGAR = r_cigar_5p_s10, read_length = calculate_read_length(r_cigar_5p_s10)
  )
  read_stats_3p_s10 <- list(
    POS = 19, CIGAR = r_cigar_3p_s10, read_length = calculate_read_length(r_cigar_3p_s10)
  )
  expect_equal(get_fragment_size(read_stats_5p_s10, read_stats_3p_s10), 25)

  # Soft clipping 11
  r_cigar_5p_s11 <- "4M5S"
  r_cigar_3p_s11 <- "4S5M"
  read_stats_5p_s11 <- list(
    POS = 1, CIGAR = r_cigar_5p_s11, read_length = calculate_read_length(r_cigar_5p_s11)
  )
  read_stats_3p_s11 <- list(
    POS = 5, CIGAR = r_cigar_3p_s11, read_length = calculate_read_length(r_cigar_3p_s11)
  )
  expect_equal(get_fragment_size(read_stats_5p_s11, read_stats_3p_s11), 9)

  # Soft clipping 12
  r_cigar_5p_s12 <- "4M2S"
  r_cigar_3p_s12 <- "5S4M"
  read_stats_5p_s12 <- list(
    POS = 1, CIGAR = r_cigar_5p_s12, read_length = calculate_read_length(r_cigar_5p_s12)
  )
  read_stats_3p_s12 <- list(
    POS = 7, CIGAR = r_cigar_3p_s12, read_length = calculate_read_length(r_cigar_3p_s12)
  )
  expect_equal(get_fragment_size(read_stats_5p_s12, read_stats_3p_s12), 10)

  # Soft clipping 13
  r_cigar_5p_s14 <- "4M5S"
  r_cigar_3p_s14 <- "2M1D7M"
  read_stats_5p_s14 <- list(
    POS = 1, CIGAR = r_cigar_5p_s14, read_length = calculate_read_length(r_cigar_5p_s14)
  )
  read_stats_3p_s14 <- list(
    POS = 5, CIGAR = r_cigar_3p_s14, read_length = calculate_read_length(r_cigar_3p_s14)
  )
  expect_equal(get_fragment_size(read_stats_5p_s14, read_stats_3p_s14), 13)

  # Soft clipping 14
  r_cigar_5p_s12 <- "4M6S"
  r_cigar_3p_s12 <- "2S2M1I5M"
  read_stats_5p_s12 <- list(
    POS = 1, CIGAR = r_cigar_5p_s12, read_length = calculate_read_length(r_cigar_5p_s12)
  )
  read_stats_3p_s12 <- list(
    POS = 7, CIGAR = r_cigar_3p_s12, read_length = calculate_read_length(r_cigar_3p_s12)
  )
  expect_equal(get_fragment_size(read_stats_5p_s12, read_stats_3p_s12), 14)
})
