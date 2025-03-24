test_that("get insertion", {
  # example 1: Case where the insertion is correctly detected

  res1 <- get_insertion(
    pos      = 15,
    alt      = "ATT",
    r_pos    = 10,
    r_cigar  = "3S6M2I1M",
    r_query  = "AAAAAAAAATTA",
    r_qual   = "#########++#"
  )
  
  expect_equal(res1$base, "insertion_detected")        
  expect_equal(res1$qual, "#++")    
  
  # example 2: Case where no insertion is present in the CIGAR
  res2 <- get_insertion(
    pos      = 15,
    alt      = "ATT",
    r_pos    = 10,
    r_cigar  = "3S7M2I",
    r_query  = "AAAAAAAAAATT",
    r_qual   = "##########++"
  )
  
  expect_equal(res2$base, "no_insertion_detected")        
  expect_equal(res2$qual, "no_insertion_detected")    


  # example 3: Case where the insertion length in CIGAR does not match alt

  res3 <- get_insertion(
    pos      = 15,
    alt      = "ATT",
    r_pos    = 10,
    r_cigar  = "3S6M1I2M",
    r_query  = "AAAAAAAAATTA",
    r_qual   = "#########+##"
  )
  
  expect_equal(res3$base, "no_insertion_detected")        
  expect_equal(res3$qual, "no_insertion_detected")  


  # example 4: Case where the position does not correspond to the insertion area

  res4 <- get_insertion(
    pos      = 15,
    alt      = "ATT",
    r_pos    = 9,
    r_cigar  = "3S6M2I1M",
    r_query  = "AAAAAAAAATTA",
    r_qual   = "#########++#"
  ) 
  
  expect_equal(res4$base, NA)
  expect_equal(res4$qual, NA)
})