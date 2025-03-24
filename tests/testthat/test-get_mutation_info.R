test_that("get insertion", {
    # example 1: insertion is correctly detected
    res1 <- get_mutation_info(
        mutation_type = "insertion",
        pos      = 15,
        ref      = "A",
        alt      = "ATT",
        read_stats <- list(pos = 10, cigar = "3S6M2I1M", query = "AAAAAAAAATTA", qual = "#########++#")
    )
    
    expect_equal(res1$base, "insertion_detected")        
    expect_equal(res1$qual, "#++")   

    # example 2: deletion is correctly detected
    res2 <- get_mutation_info(
        mutation_type = "deletion",
        pos      = 9,
        ref      = "TAT",
        alt      = "T",
        read_stats <- list(pos = 1, cigar = "9M2D10M", query = "AAAAAAAAATTA", qual = "#########++#")
    )

    expect_equal(res2$base, "deletion_detected")            
    expect_equal(res2$qual, "-")  

    # example 3: SNV is correctly detected
    res3 <- get_mutation_info(
        mutation_type = "SNV",
        pos      = 50,
        ref      = "A",
        alt      = "T",
        read_stats <- list(pos = 1, cigar = "5S5D90M", query = "CTGCCGGAACATTGGGAACCACCCCACTATGTCTAAAAGTGTTTCTTTAATCCAAATCCTCAACCCCCAATTTCCCTTCCACCCGCATAAGTTGCTGTGT", 
                      qual = ",,F,,F:,,,,,,FF,,FF,:F,,,,,,,:,,:,FFF,:::FF,,,,,,F,F,F:,:,,,F,,,F,,FFF,FF,F,F,FF,:,,F,,:F,,,,,,:,,F,")
    )

    expect_setequal(res3$base, "A")
    expect_setequal(res3$qual, "F") 

    # example 4: Not the good format of mutation
    res4 <- get_mutation_info(
        mutation_type = "SNVs",
        pos      = 50,
        ref      = "A",
        alt      = "T",
        read_stats <- list(pos = 1, cigar = "5S5D90M", query = "CTGCCGGAACATTGGGAACCACCCCACTATGTCTAAAAGTGTTTCTTTAATCCAAATCCTCAACCC", 
                      qual = ",,F,,F:,,,,,,FF,,FF,:F,,,,,,,:,,:,FFF,:::FF,,,,,,F,F,F:,:,,,F,,,F,,FFF,FF,F,F,FF,:,,F,,:F,,,,,,:,,F,")
    )

    expect_setequal(res4$base, "Error: mutation_type not defined properly")
    expect_setequal(res4$qual, "Error: mutation_type not defined properly") 
}) 