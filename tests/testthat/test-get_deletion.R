test_that("get deletion", {
    # example 1 : deletion is correctly detected
    # Assumptions:
    #  - The CIGAR string contains a deletion operation (D) of length 1.
    #  - final_pos and del_rep match (e.g., "10,1"): 
    #    => final_pos = 10, del_rep = 1.
    #  - The function check_seq_rep() returns 1 (deletion found).

    chr1       <- 1
    pos1       <- 10
    ref1       <- "T"          # the deleted base in the reference
    del_info1  <- "10,1"       # final_pos = 10, del_rep = 1
    r_pos1     <- 1
    r_cigar1   <- "9M1D10M"    # The 1D should match the expected deletion
    res1       <- get_deletion(
        chr     = chr1,
        pos     = pos1,
        ref     = ref1,
        del_info= del_info1,
        r_pos   = r_pos1,
        r_cigar = r_cigar1
    )
    
    expect_equal(res1$indel, 1)               # the function found the deletion
    expect_true(is.na(res1$base))            # base remains NA
    expect_true(is.na(res1$qual))            # qual remains NA
    
    
    # example 2: Case where final_pos is not numeric
    # Assumptions:
    #  - del_info = "NaN,1" => final_pos = NA
    #    => Expectation: c_indel = "no_del_found_in_ref_genome".

    chr2       <- 1
    pos2       <- 10
    ref2       <- "T"
    del_info2  <- "NaN,1"      # final_pos = NA => invalid case
    r_pos2     <- 1
    r_cigar2   <- "9M1D10M"    # Same CIGAR, irrelevant as the position is incorrect
    res2       <- get_deletion(
        chr     = chr2,
        pos     = pos2,
        ref     = ref2,
        del_info= del_info2,
        r_pos   = r_pos2,
        r_cigar = r_cigar2
    )
    
    expect_equal(res2$indel, "no_del_found_in_ref_genome")
    expect_true(is.na(res2$base))
    expect_true(is.na(res2$qual))
    
    
    # example 3: Case where there is no deletion in the CIGAR
    # Assumptions:
    #  - CIGAR = "20M" => no "D"
    #  - del_info = "10,1" => a deletion of length 1 is expected
    #    => Expectation: "No_support_del" (default c_indel value).

    chr3       <- 1
    pos3       <- 10
    ref3       <- "A"
    del_info3  <- "10,1"
    r_pos3     <- 1
    r_cigar3   <- "20M"       # No D => no deletion supported
    res3       <- get_deletion(
        chr     = chr3,
        pos     = pos3,
        ref     = ref3,
        del_info= del_info3,
        r_pos   = r_pos3,
        r_cigar = r_cigar3
    )
    
    expect_equal(res3$indel, "No_support_del")
    expect_true(is.na(res3$base))
    expect_true(is.na(res3$qual))
    
    
    # example 4: Case where the deletion length does not match the expected length
    # Assumptions:
    #  - CIGAR = "9M2D10M" => deletion of length 2
    #  - del_info = "10,1" => deletion of length 1 is expected
    #    => Expectation: "No_support_del" (length mismatch).

    chr4       <- 1
    pos4       <- 10
    ref4       <- "T"         # one base in the reference
    del_info4  <- "10,1"      # expecting a deletion of length 1
    r_pos4     <- 1
    r_cigar4   <- "9M2D10M"   # deletion is of length 2, not 1
    res4       <- get_deletion(
        chr     = chr4,
        pos     = pos4,
        ref     = ref4,
        del_info= del_info4,
        r_pos   = r_pos4,
        r_cigar = r_cigar4
    )
    
    expect_equal(res4$indel, "No_support_del")
    expect_true(is.na(res4$base))
    expect_true(is.na(res4$qual))
})