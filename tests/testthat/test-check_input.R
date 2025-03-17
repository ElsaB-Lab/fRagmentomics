test_that("check_input correctly handles existing and missing files", {
    # ---------------------------------------------------------------------------
    # Create temporary files
    # ---------------------------------------------------------------------------
    test_dir <- file.path(tempdir(), "test_data_check_input")
    if (!dir.exists(test_dir)) dir.create(test_dir, recursive = TRUE)
    
    # Link to the different files
    mut_missing       <- file.path(test_dir, "missing.vcf")
    mut_existing      <- file.path(test_dir, "existing.vcf.gz")
    bam_existing      <- file.path(test_dir, "existing.bam")
    bam_missing       <- file.path(test_dir, "missing.bam")
    bam_no_index      <- file.path(test_dir, "bam_no_index.bam")
    fasta_missing     <- file.path(test_dir, "missing.fa")
    fasta_no_index    <- file.path(test_dir, "no_index.fa")
    fasta_with_index  <- file.path(test_dir, "with_index.fa")

    # ---------------------------------------------------------------------------
    # Prepare fake files
    # ---------------------------------------------------------------------------
    # Create 2 fastas
    writeLines(">chr1\nAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT", fasta_no_index)
    writeLines(">chr1\nAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT", fasta_with_index)

    # Index fasta_with_index
    Rsamtools::indexFa(fasta_with_index)

    # Creat a sam file
    sam_file <- file.path(tempdir(), "test.sam")
    writeLines(c(
        "@HD\tVN:1.0\tSO:unsorted",
        "@SQ\tSN:chr1\tLN:1000",
        "read1\t0\tchr1\t1\t255\t10M\t*\t0\t0\tAGCTAGCTAG\t*"
    ), sam_file)

    # Conversion from sam to bam 
    sorted_bam <- file.path(test_dir, "sorted")

    asBam(sam_file, sorted_bam)

    # Generate existing BAM and without index
    bam_existing <- paste0(sorted_bam, ".bam")
    bam_no_index <- file.path(test_dir, "bam_no_index.bam")
    file.copy(bam_existing, bam_no_index)

    # Create index for the bam existing
    indexBam(bam_existing)

    # Create BAMs, BAI and VCF
    file.create(mut_existing)

    # ---------------------------------------------------------------------------
    # Test check_input
    # ---------------------------------------------------------------------------
    # Existing mut file
    expect_silent(
        check_input(
            mut_existing,           
            bam_existing,
            fasta_with_index
        )
    )
    
    # Missing mut file
    expect_error(
        check_input(
            mut_missing,           
            bam_existing,
            fasta_with_index
        ),
    "Error: The Mutation file does not exist",
    fixed = TRUE
    )

    # Missing BAM
    expect_error(
        check_input(
            mut_existing,
            bam_missing,          
            fasta_with_index
        ),
        "Error: The BAM file does not exist",
        fixed = TRUE
    )

    # Missing Fasta
    expect_error(
        check_input(
            mut_existing,
            bam_existing,          
            fasta_missing
        ),
        "Error: The FASTA file does not exist",
        fixed = TRUE
    )

    # Missing BAI
    bai_path <- paste0(bam_no_index, ".bai")
    if (file.exists(bai_path)) file.remove(bai_path)
    expect_message(
        check_input(
            mut_existing,
            bam_no_index,
            fasta_with_index
        ),
        "Creating BAM index..."
    )

    # Missing Fasta index 
    fai_path <- paste0(fasta_no_index, ".fai")
    if (file.exists(fai_path)) file.remove(fai_path)
    expect_message(
        check_input(
            mut_existing,
            bam_existing,
            fasta_no_index
        ),
        "Creating FASTA index..."
    )
    # ---------------------------------------------------------------------------
    # Clean test folder
    # ---------------------------------------------------------------------------
    unlink(test_dir, recursive = TRUE)
})
