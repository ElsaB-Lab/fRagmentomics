# Helper: minimal two-row SAM-like data.frame for a fragment
.make_df_sam <- function(fragment = "fragA",
                         flag1 = 0L, flag2 = 16L,
                         tlen1 = -175L, tlen2 = -175L,
                         cigar1 = "10M", cigar2 = "10M",
                         pos1 = 100L, pos2 = 200L,
                         seq1 = "ACGTACGTAA", seq2 = "TTTTCCCCGG",
                         qual1 = "FFFFFFFFFF", qual2 = "FFFFFFFFFF") {
    data.frame(
        QNAME = c(fragment, fragment), # used by the function as first column
        FLAG = c(flag1, flag2),
        MAPQ = c(60L, 60L),
        TLEN = c(tlen1, tlen2),
        CIGAR = c(cigar1, cigar2),
        POS = c(pos1, pos2),
        SEQ = c(seq1, seq2),
        QUAL = c(qual1, qual2),
        stringsAsFactors = FALSE
    )
}

# ---- extract_fragment_features() tests ---------------------------------------

test_that("Full feature extraction with TLEN, softclip counts and 5p/3p bases", {
    df <- .make_df_sam()

    # Mock all collaborators to keep this a pure orchestration test
    res <- with_mocked_bindings(
        {
            fRagmentomics:::extract_fragment_features(
                df_sam = df,
                fragment_name = "fragA",
                sample_id = "S1",
                chr = "chr1", pos = 12345L, ref = "A", alt = "T",
                report_tlen = TRUE,
                report_softclip = TRUE,
                report_5p_3p_bases_fragment = 2L,
                remove_softclip = FALSE,
                fasta_fafile = NULL, fasta_seq = NULL,
                input_mutation_info = "chr1:12345:A>T"
            )
        },
        # 1) QC passes
        process_fragment_reads_qc = function(df_fragment_reads, chr) "",
        # 2) Define one forward (5p) and one reverse (3p); row1 = forward, row2 = reverse
        bamFlagAsBitMatrix = function(flag_vec) {
            # two rows, columns as used by the function
            cbind(
                isFirstMateRead = c(TRUE, FALSE),
                isMinusStrand   = c(FALSE, TRUE)
            )
        },
        # 3) Base/qual/mutation per read
        get_base_basq_mstat_from_read = function(chr, pos, ref, alt, read_stats, fasta_fafile, fasta_seq) {
            list(base = "A", basq = "30", mstat = "MUT") # constant for simplicity
        },
        # 4) Fragment size
        get_fragment_size = function(read_stats_5p, read_stats_3p) 150L,
        # 5) Fragment status aggregation
        get_mutation_status_of_fragment = function(mstat_5p, mstat_3p) {
            list(Simple = "DISCORDANT", Detail = "MUT|MUT")
        },
        # 6) Optional 5p/3p base strings
        get_fragment_bases_5p_3p = function(n, seq5, seq3, q5, q3) {
            list(
                fragment_bases_5p = "AA",
                fragment_bases_3p = "TT",
                fragment_basqs_5p = "40 41",
                fragment_basqs_3p = "42 43"
            )
        },
        # 7) Optional softclip counts
        get_fragment_bases_5p_3p_softclip = function(cigar5, cigar3) {
            list(nb_softclip_5p = 2L, nb_softclip_3p = 3L)
        }
    )

    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 1L)
    expect_true("Sample_Id" %in% names(res))
    expect_equal(res$Sample_Id, "S1")
    expect_equal(res$Fragment_QC, "OK")
    expect_equal(res$Chromosome, "chr1")
    expect_equal(res$Position, 12345L)
    expect_equal(res$Ref, "A")
    expect_equal(res$Alt, "T")
    expect_equal(res$Input_Mutation, "chr1:12345:A>T")
    expect_equal(res$Fragment_Status_Simple, "DISCORDANT")
    expect_equal(res$Fragment_Status_Detail, "MUT|MUT")
    expect_equal(res$Fragment_Size, 150L)

    # From original rows: row1 is 5p, row2 is 3p
    expect_equal(res$FLAG_5p, df$FLAG[1])
    expect_equal(res$FLAG_3p, df$FLAG[2])
    expect_equal(res$MAPQ_5p, df$MAPQ[1])
    expect_equal(res$MAPQ_3p, df$MAPQ[2])
    expect_equal(res$CIGAR_5p, df$CIGAR[1])
    expect_equal(res$CIGAR_3p, df$CIGAR[2])
    expect_equal(res$POS_5p, df$POS[1])
    expect_equal(res$POS_3p, df$POS[2])

    # From mocks
    expect_equal(res$Read_5p_Status, "MUT")
    expect_equal(res$Read_3p_Status, "MUT")
    expect_equal(res$BASE_5p, "A")
    expect_equal(res$BASE_3p, "A")
    expect_equal(res$BASQ_5p, "30")
    expect_equal(res$BASQ_3p, "30")

    # report_tlen branch (note: code assigns numeric; do not force type)
    expect_equal(res$TLEN, abs(.make_df_sam()$TLEN[1]))

    # Optional outputs present
    expect_equal(res$Fragment_Bases_5p, "AA")
    expect_equal(res$Fragment_Bases_3p, "TT")
    expect_equal(res$Fragment_Basqs_5p, "40 41")
    expect_equal(res$Fragment_Basqs_3p, "42 43")
    expect_equal(res$Nb_Fragment_Bases_Softclip_5p, 2L)
    expect_equal(res$Nb_Fragment_Bases_Softclip_3p, 3L)

    # VAF always added, NA_real_
    expect_true(is.na(res$VAF))
    expect_type(res$VAF, "double")
})

test_that("QC failure returns the placeholder row from create_empty_fragment_row", {
    df <- .make_df_sam()
    placeholder <- data.frame(
        Chromosome = "chrX",
        Position = 999L,
        Ref = "G", Alt = "C",
        Input_Mutation = "chrX:999:G>C",
        Fragment_Id = "fragA",
        Fragment_QC = "QC_FAIL",
        stringsAsFactors = FALSE
    )

    res <- with_mocked_bindings(
        {
            fRagmentomics:::extract_fragment_features(
                df_sam = df,
                fragment_name = "fragA",
                sample_id = "S_FAIL",
                chr = "chrX", pos = 999L, ref = "G", alt = "C",
                report_tlen = TRUE,
                report_softclip = TRUE,
                report_5p_3p_bases_fragment = 1L,
                remove_softclip = FALSE,
                fasta_fafile = NULL, fasta_seq = NULL,
                input_mutation_info = "chrX:999:G>C"
            )
        },
        process_fragment_reads_qc = function(df_fragment_reads, chr) "QC_FAIL",
        create_empty_fragment_row = function(chr, pos, ref, alt, input_mutation_info,
                                             fragment_name, fragment_qc, sample_id,
                                             report_tlen, report_5p_3p_bases_fragment,
                                             report_softclip) {
            placeholder
        }
    )

    expect_equal(res, placeholder)
})

test_that("Stops when not exactly one first mate (invalid R1/R2 pairing)", {
    df <- .make_df_sam()

    expect_error(
        with_mocked_bindings(
            fRagmentomics:::extract_fragment_features(
                df_sam = df,
                fragment_name = "fragA",
                sample_id = NA_character_,
                chr = "chr1", pos = 1L, ref = "A", alt = "T",
                report_tlen = FALSE, report_softclip = FALSE,
                report_5p_3p_bases_fragment = 0L, remove_softclip = FALSE,
                fasta_fafile = NULL, fasta_seq = NULL,
                input_mutation_info = "chr1:1:A>T"
            ),
            process_fragment_reads_qc = function(df_fragment_reads, chr) "",
            bamFlagAsBitMatrix = function(flag_vec) {
                cbind(
                    isFirstMateRead = c(TRUE, TRUE), # <- invalid (sum != 1)
                    isMinusStrand   = c(FALSE, TRUE)
                )
            }
        ),
        "not a valid R1/R2 pair"
    )
})

test_that("Stops when not exactly one forward and one reverse read", {
    df <- .make_df_sam()

    expect_error(
        with_mocked_bindings(
            fRagmentomics:::extract_fragment_features(
                df_sam = df,
                fragment_name = "fragA",
                sample_id = NA_character_,
                chr = "chr1", pos = 1L, ref = "A", alt = "T",
                report_tlen = FALSE, report_softclip = FALSE,
                report_5p_3p_bases_fragment = 0L, remove_softclip = FALSE,
                fasta_fafile = NULL, fasta_seq = NULL,
                input_mutation_info = "chr1:1:A>T"
            ),
            process_fragment_reads_qc = function(df_fragment_reads, chr) "",
            bamFlagAsBitMatrix = function(flag_vec) {
                cbind(
                    isFirstMateRead = c(TRUE, FALSE),
                    isMinusStrand   = c(FALSE, FALSE) # <- invalid (sum != 1)
                )
            }
        ),
        "does not have one forward and one reverse read"
    )
})

test_that("remove_softclip=TRUE: early return when 5p becomes empty", {
    # use sentinel CIGAR to branch inside our mock
    df <- .make_df_sam(cigar1 = "EMPTY5", cigar2 = "OK")

    placeholder <- data.frame(
        Chromosome = "chr1", Position = 123L, Ref = "A", Alt = "T",
        Input_Mutation = "chr1:123:A>T",
        Fragment_Id = "fragA",
        Fragment_QC = "Read empty after softclip trim",
        stringsAsFactors = FALSE
    )

    res <- with_mocked_bindings(
        {
            fRagmentomics:::extract_fragment_features(
                df_sam = df,
                fragment_name = "fragA",
                sample_id = "Sx",
                chr = "chr1", pos = 123L, ref = "A", alt = "T",
                report_tlen = TRUE, report_softclip = TRUE,
                report_5p_3p_bases_fragment = 2L, remove_softclip = TRUE,
                fasta_fafile = NULL, fasta_seq = NULL,
                input_mutation_info = "chr1:123:A>T"
            )
        },
        process_fragment_reads_qc = function(df_fragment_reads, chr) "",
        bamFlagAsBitMatrix = function(flag_vec) {
            cbind(isFirstMateRead = c(TRUE, FALSE), isMinusStrand = c(FALSE, TRUE))
        },
        remove_softclip = function(read_stats) {
            if (identical(read_stats$CIGAR, "EMPTY5")) {
                # 5p becomes empty -> triggers early return
                list(SEQ = "", QUAL = "", CIGAR = "", read_length = 0L)
            } else {
                # pass-through for 3p
                list(
                    SEQ = read_stats$SEQ, QUAL = read_stats$QUAL,
                    CIGAR = read_stats$CIGAR, read_length = nchar(read_stats$SEQ)
                )
            }
        },
        create_empty_fragment_row = function(chr, pos, ref, alt, input_mutation_info,
                                             fragment_name, fragment_qc, sample_id,
                                             report_tlen, report_5p_3p_bases_fragment,
                                             report_softclip) {
            placeholder
        }
    )

    expect_equal(res, placeholder)
})

test_that("remove_softclip=TRUE: early return when 3p becomes empty", {
    df <- .make_df_sam(cigar1 = "OK", cigar2 = "EMPTY3")

    placeholder <- data.frame(
        Chromosome = "chr1", Position = 456L, Ref = "G", Alt = "C",
        Input_Mutation = "chr1:456:G>C",
        Fragment_Id = "fragA",
        Fragment_QC = "Read empty after softclip trim",
        stringsAsFactors = FALSE
    )

    res <- with_mocked_bindings(
        {
            fRagmentomics:::extract_fragment_features(
                df_sam = df,
                fragment_name = "fragA",
                sample_id = "Sy",
                chr = "chr1", pos = 456L, ref = "G", alt = "C",
                report_tlen = TRUE, report_softclip = TRUE,
                report_5p_3p_bases_fragment = 1L, remove_softclip = TRUE,
                fasta_fafile = NULL, fasta_seq = NULL,
                input_mutation_info = "chr1:456:G>C"
            )
        },
        process_fragment_reads_qc = function(df_fragment_reads, chr) "",
        bamFlagAsBitMatrix = function(flag_vec) {
            cbind(isFirstMateRead = c(TRUE, FALSE), isMinusStrand = c(FALSE, TRUE))
        },
        remove_softclip = function(read_stats) {
            if (identical(read_stats$CIGAR, "EMPTY3")) {
                list(SEQ = "", QUAL = "", CIGAR = "", read_length = 0L)
            } else {
                list(
                    SEQ = read_stats$SEQ, QUAL = read_stats$QUAL,
                    CIGAR = read_stats$CIGAR, read_length = nchar(read_stats$SEQ)
                )
            }
        },
        create_empty_fragment_row = function(chr, pos, ref, alt, input_mutation_info,
                                             fragment_name, fragment_qc, sample_id,
                                             report_tlen, report_5p_3p_bases_fragment,
                                             report_softclip) {
            placeholder
        }
    )

    expect_equal(res, placeholder)
})

test_that("TLEN reporting: warning when TLEN is missing; no Sample_Id column when NA", {
    df <- .make_df_sam(tlen1 = NA_integer_, tlen2 = NA_integer_)

    # Prepare minimal mocks to pass through main pipeline,
    # with report_5p_3p_bases_fragment = 0 and report_softclip = FALSE
    expect_warning(
        res <- with_mocked_bindings(
            {
                fRagmentomics:::extract_fragment_features(
                    df_sam = df,
                    fragment_name = "fragA",
                    sample_id = NA_character_,
                    chr = "chr2", pos = 888L, ref = "C", alt = "A",
                    report_tlen = TRUE,
                    report_softclip = FALSE,
                    report_5p_3p_bases_fragment = 0L,
                    remove_softclip = FALSE,
                    fasta_fafile = NULL, fasta_seq = NULL,
                    input_mutation_info = "chr2:888:C>A"
                )
            },
            process_fragment_reads_qc = function(df_fragment_reads, chr) "",
            bamFlagAsBitMatrix = function(flag_vec) {
                cbind(isFirstMateRead = c(TRUE, FALSE), isMinusStrand = c(FALSE, TRUE))
            },
            get_base_basq_mstat_from_read = function(chr, pos, ref, alt, read_stats, fasta_fafile, fasta_seq) {
                list(base = "C", basq = "25", mstat = "REF")
            },
            get_fragment_size = function(read_stats_5p, read_stats_3p) 101L,
            get_mutation_status_of_fragment = function(mstat_5p, mstat_3p) {
                list(Simple = "REF", Detail = "REF|REF")
            }
        ),
        "TLEN is NULL or missing"
    )

    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 1L)
    expect_false("Sample_Id" %in% names(res)) # sample_id = NA -> column omitted
    expect_true("TLEN" %in% names(res))
    expect_equal(res$TLEN, "Warning: TLEN is NULL") # as specified by implementation

    # Since report_5p_3p_bases_fragment = 0 and report_softclip = FALSE,
    # no extra columns should be present
    expect_false(any(c(
        "Fragment_Bases_5p", "Fragment_Bases_3p",
        "Fragment_Basqs_5p", "Fragment_Basqs_3p"
    ) %in% names(res)))
    expect_false(any(c(
        "Nb_Fragment_Bases_Softclip_5p",
        "Nb_Fragment_Bases_Softclip_3p"
    ) %in% names(res)))
})
