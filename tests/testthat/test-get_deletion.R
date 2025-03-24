test_that("get deletion", {
    # example 1 : deletion is correctly detected
    res1       <- get_deletion(
        pos     = 9,                       # Position target the nucleotide before the deletion
        ref     = "TAT",
        r_pos   = 1,
        r_cigar = "9M2D10M" 
    )
    
    expect_equal(res1$base, "deletion_detected")            
    expect_equal(res1$qual, "-")          
    
    # example 2: Case where the read doesn't carry the deletion
    res2       <- get_deletion(
        pos     = 9,
        ref     = "TA",
        r_pos   = 1,
        r_cigar = "15M1D4M"
    )
    
    expect_equal(res2$base, "no_deletion_detected")
    expect_equal(res2$qual, "no_deletion_detected")
    
    # example 3: Case where the deletion length does not match the expected length
    res3       <- get_deletion(
        pos     = 9,
        ref     = "TAT",
        r_pos   = 1,
        r_cigar = "8M2D"
    )
    
    expect_equal(res3$base, "no_deletion_detected")
    expect_equal(res3$qual, "no_deletion_detected")

    # example 4: Case where the read doesn't cover the position of the deletion
    res3       <- get_deletion(
        pos     = 9,
        ref     = "TAT",
        r_pos   = 1,
        r_cigar = "9M"
    )
    
    expect_equal(res3$base, NA)
    expect_equal(res3$qual, NA)
})