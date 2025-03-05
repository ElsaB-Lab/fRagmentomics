test_that("check_input correctly handles existing and missing files", {
    # Create a temporary test directory
    test_dir <- file.path(tempdir(), "test_data")
    if (!dir.exists(test_dir)) dir.create(test_dir, recursive = TRUE)

    bam_existing <- file.path(test_dir, "existing.bam")
    bam_missing <- file.path(test_dir, "missing.bam")
    fasta_missing <- file.path(test_dir, "missing.fa")
    fasta_no_index <- file.path(test_dir, "no_index.fa")
    fasta_fai_existing <- file.path(test_dir, "all_good.fa")

    # Write valid FASTA files
    writeLines(">chr1\nAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT", fasta_no_index)
    writeLines(">chr1\nAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT", fasta_fai_existing)

    # Create a dummy BAM file
    file.create(bam_existing)

    # Ensure a valid FASTA index exists for `all_good.fa`
    Rsamtools::indexFa(fasta_fai_existing)

    # Case 1: Fasta, BAM, and Index existing -> Expect silent execution
    expect_silent(check_input(bam_existing, fasta_fai_existing))

    # Case 2: BAM file is missing -> Expect error
    expect_error(check_input(bam_missing, fasta_fai_existing), "Error: The BAM file does not exist")

    # Case 3: FASTA file is missing -> Expect error
    expect_error(check_input(bam_existing, fasta_missing), "Error: The FASTA file does not exist")

    # Case 4: FASTA file exists but no index
    # Ensure the FASTA index does not exist before testing
    fai_path <- paste0(fasta_no_index, ".fai")
    if (file.exists(fai_path)) file.remove(fai_path)
    expect_message(check_input(bam_existing, fasta_no_index))

    # Cleanup test files
    unlink(test_dir, recursive = TRUE)
})
