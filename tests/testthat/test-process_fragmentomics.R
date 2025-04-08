test_that("Process fRagmentomics", {
    mut_test <- system.file("testdata/mutations/", "mutations_cfdna-test-01_chr17_7576000_7579000.tsv", package = "fRagmentomics")
    bam_test <- system.file("testdata/bam/", "cfdna-test-01_chr17_7576000_7579000.bam", package = "fRagmentomics")
    fasta_test <- system.file("testdata/fasta/hg19/", "hg19_chr17_7576000_7579000.fa", package = "fRagmentomics")

    # Read the mutation file
    mutation_file_expected <- read.table(
        mut_test,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
    )

    # Execute the function with sample_id = NA
    res1 <- process_fragmentomics(
        mut = mut_test,
        bam = bam_test,
        fasta = fasta_test,
        sample_id = NA,
        neg_offset_mate_search = -200,
        pos_offset_mate_search = 200,
        one_based = TRUE,
        flag_keep = 0x03,
        flag_remove = 0x900,
        report_tlen = TRUE,
        report_softclip = TRUE,
        report_5p_3p_bases_fragment = 5,
        tmp_folder = tempdir(),
        output_folder = ".", # file will be created in the current directory
        n_cores = 1
    )

    # Check that each row of res1 matches a line in the mutation file
    for (i in seq_len(nrow(res1))) {
        expect_true(
            any(
                mutation_file_expected$CHROM == res1$Chromosome[i] &
                    mutation_file_expected$POS == res1$Position[i] &
                    mutation_file_expected$REF == res1$Ref[i] &
                    mutation_file_expected$ALT == res1$Alt[i]
            )
        )
    }

    # Check that we have the expected columns
    expected_cols <- c(
        "Chromosome", "Position", "Ref", "Alt",
        "Fragment_QC", "Fragment_Mutated", "Absolute_size",
        "Inner_distance", "MAPQ_5p", "MAPQ_3p",
        "BASE_5p", "BASE_3p", "BASQ_5p", "BASQ_3p"
    )

    expect_true(
        all(expected_cols %in% colnames(res1)),
        info = "Some expected columns are missing in the output dataframe."
    )

    # Now delete the output file created
    created_files <- list.files(
        path = ".",
        pattern = "^[0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]{2}:[0-9]{2}:[0-9]{2}_fRagmentomics_output\\.tsv$"
    )

    # We expect at least one file with that naming; fail if not found
    expect_true(
        length(created_files) > 0,
        info = "No output file found matching the pattern 'YYYY-mm-dd_HH:MM:SS_fRagmentomics_output.tsv'."
    )

    # Remove each matching file and verify it was removed
    for (f in created_files) {
        file_removed <- file.remove(f)
        expect_true(
            file_removed,
            info = paste("Failed to remove output file", f)
        )
        expect_false(
            file.exists(f),
            info = paste("File still exists after removal attempt:", f)
        )
    }
})
