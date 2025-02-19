test_that("get pos indel from read", {
  # example 1 : sample ORD-1014010-01 and fragment 07102AuSP:0602
  r_pos1 <- 7578167
  r_cigar1 <- "125M1D14M5S"
  r_pos2 <- 7578196
  r_cigar2 <- "17S96M1D31M"
  
  d_pos_i_pos1 <- get_pos_indel_from_read(r_pos=r_pos1, r_cigar=r_cigar1) 
  d_pos_i_pos2 <- get_pos_indel_from_read(r_pos=r_pos2, r_cigar=r_cigar2)

  expect_true(setequal(d_pos_i_pos1$d_pos, 7578292))
  expect_true(setequal(d_pos_i_pos1$i_pos, NULL))
  expect_true(setequal(d_pos_i_pos2$d_pos, 7578292))
  expect_true(setequal(d_pos_i_pos2$i_pos, NULL))

  # example 2 :
  r_pos3 <- 7578027
  r_cigar3 <- "144M" 
  r_pos4 <- 7578078
  r_cigar4 <- "119M1D25M"
 
  d_pos_i_pos3 <- get_pos_indel_from_read(r_pos=r_pos3, r_cigar=r_cigar3)
  d_pos_i_pos4 <- get_pos_indel_from_read(r_pos=r_pos4, r_cigar=r_cigar4)

  expect_true(setequal(d_pos_i_pos3$d_pos, NULL))
  expect_true(setequal(d_pos_i_pos3$i_pos, NULL))
  expect_true(setequal(d_pos_i_pos4$d_pos, 7578197))
  expect_true(setequal(d_pos_i_pos4$i_pos, NULL))


  # example 3 : 
  r_pos5 <- 1
  r_cigar5 <- "5S5D10M10I70M"
  
  d_pos_i_pos5 <- get_pos_indel_from_read(r_pos=r_pos5, r_cigar=r_cigar5)
  
  expect_true(setequal(d_pos_i_pos5$d_pos, c(1,2,3,4,5)))
  expect_true(setequal(d_pos_i_pos5$i_pos, c(16,16,16,16,16,16,16,16,16,16)))


})
