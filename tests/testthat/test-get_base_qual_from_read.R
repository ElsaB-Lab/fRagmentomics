test_that("get base and qual from read", {
  # example 1
  r_pos1 <- 1
  pos1 <- 3
  r_cigar1 <- "2M1D50M"
  r_query1 <- "CTGCCGGAACATTGGGAACCACCCCACTATGTCTAAAAGTGTTTCTTTAATC"
  r_qual1 <- ",,F,,F:,,,,,,FF,,FF,:F,,,,,,,:,,:,FFF,:::FF,,,,,,F,F"
  base_qual1 <- get_base_qual_from_read(
    pos = pos1,
    r_pos = r_pos1,
    r_cigar = r_cigar1,
    r_query = r_query1,
    r_qual = r_qual1
  )

  expect_setequal(base_qual1$base, "-")
  expect_setequal(base_qual1$qual, NA)

  # example 2
  r_pos2 <- 1
  pos2 <- 100
  r_cigar2 <- "2M1D50M"
  r_query2 <- "CTGCCGGAACATTGGGAACCACCCCACTATGTCTAAAAGTGTTTCTTTAATC"
  r_qual2 <- ",,F,,F:,,,,,,FF,,FF,:F,,,,,,,:,,:,FFF,:::FF,,,,,,F,F"
  base_qual2 <- get_base_qual_from_read(
    pos = pos2,
    r_pos = r_pos2,
    r_cigar = r_cigar2,
    r_query = r_query2,
    r_qual = r_qual2
  )

  expect_setequal(base_qual2$base, NA)
  expect_setequal(base_qual2$qual, NA)

  # example 3
  r_pos3 <- 1
  pos3 <- 50
  r_cigar3 <- "5S95M"
  r_query3 <- "CTGCCGGAACATTGGGAACCACCCCACTATGTCTAAAAGTGTTTCTTTAATCCAAATCCTCAACCCCCAATTTCCCTTCCACCCGCATAAGTTGCTGTGT"
  r_qual3 <- ",,F,,F:,,,,,,FF,,FF,:F,,,,,,,:,,:,FFF,:::FF,,,,,,F,F,F:,:,,,F,,,F,,FFF,FF,F,F,FF,:,,F,,:F,,,,,,:,,F,"
  base_qual3 <- get_base_qual_from_read(pos = pos3, r_pos = r_pos3, r_cigar = r_cigar3, r_query = r_query3, r_qual = r_qual3)

  expect_setequal(base_qual3$base, "A")
  expect_setequal(base_qual3$qual, ":")

  # example 4
  r_pos4 <- 1
  pos4 <- 50
  r_cigar4 <- "5S5D90M"
  r_query4 <- "CTGCCGGAACATTGGGAACCACCCCACTATGTCTAAAAGTGTTTCTTTAATCCAAATCCTCAACCCCCAATTTCCCTTCCACCCGCATAAGTTGCTGTGT"
  r_qual4 <- ",,F,,F:,,,,,,FF,,FF,:F,,,,,,,:,,:,FFF,:::FF,,,,,,F,F,F:,:,,,F,,,F,,FFF,FF,F,F,FF,:,,F,,:F,,,,,,:,,F,"
  base_qual4 <- get_base_qual_from_read(pos = pos4, r_pos = r_pos4, r_cigar = r_cigar4, r_query = r_query4, r_qual = r_qual4)

  expect_setequal(base_qual4$base, "A")
  expect_setequal(base_qual4$qual, "F")

  # example 5
  r_pos5 <- 1
  pos5 <- 50
  r_cigar5 <- "5S5D10M10I70M"
  r_query5 <- "CTGCCGGAACATTGGGAACCACCCCACTATGTCTAAAAGTGTTTCTTTAATCCAAATCCTCAACCCCCAATTTCCCTTCCACCCGCATAAGTTGC"
  r_qual5 <- ",,F,,F:,,,,,,FF,,FF,:F,,,,,,,:,,:,FFF,:::FF,,,,,,F,F,F:,:,,,F,,,F,,FFF,FF,F,F,FF,:,,F,,:F,,,,,,"
  base_qual5 <- get_base_qual_from_read(pos = pos5, r_pos = r_pos5, r_cigar = r_cigar5, r_query = r_query5, r_qual = r_qual5)

  expect_setequal(base_qual5$base, "T")
  expect_setequal(base_qual5$qual, ",")

  # example 6
  r_pos6 <- 1
  pos6 <- 50
  r_cigar6 <- "60M10D10I20S"
  r_query6 <- "CTGCCGGAACATTGGGAACCACCCCACTATGTCTAAAAGTGTTTCTTTAATCCAAATCCTCAACCCCCAATTTCCCTTCCACCCGCATAA"
  r_qual6 <- ",,F,,F:,,,,,,FF,,FF,:F,,,,,,,:,,:,FFF,:::FF,,,,,,F,F,F:,:,,,F,,,F,,FFF,FF,F,F,FF,:,,F,,:F,"
  base_qual6 <- get_base_qual_from_read(pos = pos6, r_pos = r_pos6, r_cigar = r_cigar6, r_query = r_query6, r_qual = r_qual6)

  expect_setequal(base_qual6$base, "A")
  expect_setequal(base_qual6$qual, "F")

  # example 7
  r_pos7 <- 1
  pos7 <- 50
  r_cigar7 <- "49M1I50M"
  r_query7 <- "CTGCCGGAACATTGGGAACCACCCCACTATGTCTAAAAGTGTTTCTTTAATCCAAATCCTCAACCCCCAATTTCCCTTCCACCCGCATAAGTTGCTGTGT"
  r_qual7 <- ",,F,,F:,,,,,,FF,,FF,:F,,,,,,,:,,:,FFF,:::FF,,,,,,F!F,F:,:,,,F,,,F,,FFF,FF,F,F,FF,:,,F,,:F,,,,,,:,,F,"
  base_qual7 <- get_base_qual_from_read(pos = pos7, r_pos = r_pos7, r_cigar = r_cigar7, r_query = r_query7, r_qual = r_qual7)

  expect_setequal(base_qual7$base, "T")
  expect_setequal(base_qual7$qual, "!")

  # example 8
  r_pos8 <- 1
  pos8 <- 105
  r_cigar8 <- "10S80M10S"
  r_query8 <- "CTGCCGGAACATTGGGAACCACCCCACTATGTCTAAAAGTGTTTCTTTAATCCAAATCCTCAACCCCCAATTTCCCTTCCACCCGCATAAGTTGCTGTGT"
  r_qual8 <- ",,F,,F:,,,,,,FF,,FF,:F,,,,,,,:,,:,FFF,:::FF,,,,,,F,F,F:,:,,,F,,,F,,FFF,FF,F,F,FF,:,,F,,:F,,,,,,:,,F,"
  base_qual8 <- get_base_qual_from_read(pos = pos8, r_pos = r_pos8, r_cigar = r_cigar8, r_query = r_query8, r_qual = r_qual8)

  expect_setequal(base_qual8$base, NA)
  expect_setequal(base_qual8$qual, NA)


  # example 9
  r_pos9 <- 1
  pos9 <- 50
  r_cigar9 <- "5S5I10X90="
  r_query9 <- "CTGCCGGAACATTGGGAACCACCCCACTATGTCTAAAAGTGTTTCTTTAATCCAAATCCTCAACCCCCAATTTCCCTTCCACCCGCATAAGTTGCTGTGTTGGGGCCAGA"
  r_qual9 <- ",,F,,F:,,,,,,FF,,FF,:F,,,,,,,:,,:,FFF,:::FF,,,,,,F,F,F:,:,,,F,,,F,,FFF,FF,F,F,FF,:,,F,,:F,,,,,,:,,F,,FFFFFF,,,"
  base_qual9 <- get_base_qual_from_read(pos = pos9, r_pos = r_pos9, r_cigar = r_cigar9, r_query = r_query9, r_qual = r_qual9)

  expect_setequal(base_qual9$base, "T")
  expect_setequal(base_qual9$qual, ",")

  # example 10 -> Issue when the mutation is the first base of the soft clipping
  r_pos10 <- 6
  pos10 <- 9
  r_cigar10 <- "5S3M3S"
  r_query10 <- "ATCGAATCGTX"
  r_qual10 <- "FFFFFFFF!FX"
  base_qual10 <- get_base_qual_from_read(pos = pos10, r_pos = r_pos10, r_cigar = r_cigar10, r_query = r_query10, r_qual = r_qual10)

  expect_setequal(base_qual10$base, NA)
  expect_setequal(base_qual10$qual, NA)

  # example 10 -> Test if SNV on the last position before softclipping
  r_pos11 <- 6
  pos11 <- 8
  r_cigar11 <- "5S3M3S"
  r_query11 <- "ATCGAATCGTX"
  r_qual11 <- "FFFFFFF!FFX"
  base_qual11 <- get_base_qual_from_read(pos = pos11, r_pos = r_pos11, r_cigar = r_cigar11, r_query = r_query11, r_qual = r_qual11)

  expect_setequal(base_qual11$base, "C")
  expect_setequal(base_qual11$qual, "!")
})
