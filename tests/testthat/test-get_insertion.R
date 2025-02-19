test_that("get insertion", {
  
  # example 1: Case where the insertion is correctly detected
  # Assumptions:
  #  - The CIGAR contains an I operation of length 3.
  #  - alt = "ATG" (length = 3)
  #  - Internal functions (define_ins_rep, check_seq_rep, find_c_qual) return valid results.

  pos1       <- 10
  alt1       <- "ATG"        # the expected insertion
  r_pos1     <- 1
  r_cigar1   <- "9M3I10M"    # CIGAR with 3 bases inserted after 9 matches
  r_query1   <- "AAAAAAAAAATGTTTTTTTTTT"
  r_qual1    <- "#########+++##########"
  res1 <- get_insertion(
    pos      = pos1,
    alt      = alt1,
    r_pos    = r_pos1,
    r_cigar  = r_cigar1,
    r_query  = r_query1,
    r_qual   = r_qual1
  )
  
  expect_equal(res1$indel, 1)          # the function detected the insertion
  expect_true(is.na(res1$base))        # base remains NA
  expect_equal(nchar(res1$qual), 3)    # quality matches insertion length
  

  # example 2: Case where no insertion is present in the CIGAR
  # Assumptions:
  #  - CIGAR = "20M" => no "I"
  #  - alt = "ATG"
  #    => Expectation: "No_support_del".

  pos2       <- 10
  alt2       <- "ATG"
  r_pos2     <- 1
  r_cigar2   <- "20M"
  r_query2   <- "AAAAAAAAAAAAAAAAAAAA"
  r_qual2    <- "####################"
  res2 <- get_insertion(
    pos      = pos2,
    alt      = alt2,
    r_pos    = r_pos2,
    r_cigar  = r_cigar2,
    r_query  = r_query2,
    r_qual   = r_qual2
  )
  
  expect_equal(res2$indel, "No_support_ins")
  expect_true(is.na(res2$base))
  expect_true(is.na(res2$qual))


  # example 3: Case where the insertion length in CIGAR does not match alt
  # Assumptions:
  #  - CIGAR = "9M2I10M" => insertion of 2 bases
  #  - alt = "ATG" => 3 bases
  #    => Expectation: "No_support_ins".

  pos3       <- 10
  alt3       <- "ATG"   # 3 bases
  r_pos3     <- 1
  r_cigar3   <- "9M2I10M"
  r_query3   <- "AAAAAAAAAATTTTTTTTTT"
  r_qual3    <- "##########++++++++++"
  res3 <- get_insertion(
    pos      = pos3,
    alt      = alt3,
    r_pos    = r_pos3,
    r_cigar  = r_cigar3,
    r_query  = r_query3,
    r_qual   = r_qual3
  )
  
  expect_equal(res3$indel, "No_support_ins")
  expect_true(is.na(res3$base))
  expect_true(is.na(res3$qual))


  # example 4: Case where the position does not correspond to the insertion area
  # Assumptions:
  #  - alt_len = 3
  #  - CIGAR = "9M3I10M", insertion of 3 bases
  #  - BUT the expected position `pos4` is 50, outside the range covered by the insertion.
  #    => Expectation: "No_support_del".

  pos4       <- 50    # position too far
  alt4       <- "ATG"
  r_pos4     <- 1
  r_cigar4   <- "9M3I10M"
  r_query4   <- "AAAAAAAAAATGTTTTTTTTTT"
  r_qual4    <- "#########+++##########"
  res4 <- get_insertion(
    pos      = pos4,
    alt      = alt4,
    r_pos    = r_pos4,
    r_cigar  = r_cigar4,
    r_query  = r_query4,
    r_qual   = r_qual4
  )
  
  expect_equal(res4$indel, "No_support_ins")
  expect_true(is.na(res4$base))
  expect_true(is.na(res4$qual))
})