test_that("build fragments info table on good fragments", {
  library(stringr)
  chr <- "chr17"
  pos <- 7578190
  df_info <- build_fragments_info_table(df_sam, chr, pos)

  expect_true(setequal(colnames(df_info), c("Chromosome", "Position", "Fragment", "Fragment_Check", "Mapq_Read_1",
                                            "Mapq_Read_2", "Base_Read_1", "Base_Read_2", "Qual_Read_1", "Qual_Read_2",
                                            "TLEN", "Cigar_Read_1", "Cigar_Read_2",
					    "Read_length_1", "Read_length_2",
                                            "Startpos_read_1", "Startpos_read_2","Identity_RIGHT", "Identity_LEFT",
                                            "Bases_del", "Bases_ins", "Bases_soft_clip_left", "Bases_soft_clip_right",
                                            "Inner_distance_a", "Inner_distance_b", "Absolute_size")))
  expect_true(all(df_info$Chromosome==chr))
  expect_true(all(df_info$Position==pos))
  expect_true(setequal(df_info$Fragment, unique(df_sam$qname)))
  expect_true(all(df_info$Fragment_Check=="OK"))
})
