test_that("define_mutation_status handles different mutation types correctly", {
    
    #---------------------------------------
    # SNV: Single Nucleotide Variant (1 base change)
    #---------------------------------------
    expect_equal(define_mutation_status("A", "T"), "SNV")
    
    #---------------------------------------
    # MNP: Multiple bases changed, same length
    #---------------------------------------
    expect_equal(define_mutation_status("AG", "TC"), "MNP")
    
    #---------------------------------------
    # Deletion: REF longer than ALT
    #---------------------------------------
    expect_equal(define_mutation_status("AGT", "A"), "deletion")
    
    #---------------------------------------
    # Insertion: ALT longer than REF
    #---------------------------------------
    expect_equal(define_mutation_status("A", "ACT"), "insertion")
})
