# Helper: minimal two-row SAM-like data.frame for a fragment
.make_df_sam <- function(fragment = "fragA",
                         flag1 = 163L, flag2 = 83L,
                         tlen1 = -175L, tlen2 = -175L,
                         cigar1 = "10M", cigar2 = "10M",
                         pos1 = 100L, pos2 = 200L,
                         seq1 = "ACGTACGTAA", seq2 = "TTTTCCCCGG",
                         qual1 = "FFFFFFFFFF", qual2 = "FFFFFFFFFF",
                         rname1 = "chr1", rname2 = "chr1",
                         rnext = "=",
                         pnext1 = 200L, pnext2 = 100L) {
    data.frame(
        QNAME = c(fragment, fragment),
        FLAG = c(flag1, flag2),
        MAPQ = c(60L, 60L),
        TLEN = c(tlen1, tlen2),
        CIGAR = c(cigar1, cigar2),
        RNAME = c(rname1, rname2),
        POS = c(pos1, pos2),
        RNEXT = c(rnext, rnext),
        PNEXT = c(pnext1, pnext2),
        SEQ = c(seq1, seq2),
        QUAL = c(qual1, qual2),
        stringsAsFactors = FALSE
    )
}

# ---- extract_fragment_features() tests ---------------------------------------

test_that("Full feature extraction with TLEN, softclip counts and 5p/3p bases", {
    df <- .make_df_sam()

    res <- with_mocked_bindings(
        {
            fRagmentomics:::extract_fragment_features(
                df_sam = df,
                fragment_name = "fragA",
                sample_id = "S1",
                chr = "chr1", pos = 12345L, ref = "A", alt = "T",
                report_bam_info = TRUE,
                report_softclip = TRUE,
                report_5p_3p_bases_fragment = 2L,
                remove_softclip = FALSE, # renamed boolean flag
                fasta_fafile = NULL, fasta_seq = NULL,
                input_mutation_info = "chr1:12345:A>T"
            )
        },
        # 1) QC passes
        process_fragment_reads_qc = function(df_fragment_reads, chr) "",
        # 2) Define one forward (5p) and one reverse (3p); row1 = forward, row2 = reverse
        bamFlagAsBitMatrix = function(flag_vec) {
            cbind(
                isFirstMateRead = c(TRUE, FALSE),
                isMinusStrand   = c(FALSE, TRUE)
            )
        },
        # 3) Base/qual/mutation per read
        get_base_basq_mstat_from_read = function(chr, pos, ref, alt, read_stats, fasta_fafile, fasta_seq) {
            list(base = "A", basq = "30", mstat = "MUT")
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
        },
        # 8) Last aligned position (3p) used inside end_on_reference()
        parse_cigar = function(cigar) {
            # emulate parse_cigar output: data.frame with 'length' and 'type'
            # For "10M", the ref length contribution is 10
            data.frame(length = as.integer(10L), type = "M")
        }
    )

    # Function returns a named list (assembled by rbindlist higher up)
    expect_type(res, "list")
    expect_true(all(c(
        "Sample_Id", "Chromosome", "Position", "Ref", "Alt", "Input_Mutation", "Fragment_Id", "Fragment_QC",
        "Fragment_Status_Simple", "Fragment_Status_Detail", "Fragment_Size",
        "Read_5p_Status", "Read_3p_Status", "BASE_5p", "BASE_3p", "BASQ_5p", "BASQ_3p",
        "Position_5p", "Position_3p", "VAF",
        "POS_5p", "POS_3p", "FLAG_5p", "FLAG_3p", "MAPQ_5p", "MAPQ_3p", "CIGAR_5p", "CIGAR_3p",
        "Fragment_Bases_5p", "Fragment_Bases_3p", "Fragment_Basqs_5p", "Fragment_Basqs_3p",
        "Nb_Fragment_Bases_Softclip_5p", "Nb_Fragment_Bases_Softclip_3p", "TLEN"
    ) %in% names(res)))

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
    expect_equal(res$FLAG_5p, .make_df_sam()$FLAG[1])
    expect_equal(res$FLAG_3p, .make_df_sam()$FLAG[2])
    expect_equal(res$MAPQ_5p, .make_df_sam()$MAPQ[1])
    expect_equal(res$MAPQ_3p, .make_df_sam()$MAPQ[2])
    expect_equal(res$CIGAR_5p, .make_df_sam()$CIGAR[1])
    expect_equal(res$CIGAR_3p, .make_df_sam()$CIGAR[2])
    expect_equal(res$POS_5p, .make_df_sam()$POS[1])
    expect_equal(res$POS_3p, .make_df_sam()$POS[2])

    # From mocks
    expect_equal(res$Read_5p_Status, "MUT")
    expect_equal(res$Read_3p_Status, "MUT")
    expect_equal(res$BASE_5p, "A")
    expect_equal(res$BASE_3p, "A")
    expect_equal(res$BASQ_5p, "30")
    expect_equal(res$BASQ_3p, "30")

    # TLEN: positive absolute as integer
    expect_equal(res$TLEN, as.integer(abs(.make_df_sam()$TLEN[1])))

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

test_that("QC failure returns a single-row placeholder list with NAs and Fragment_QC message", {
    df <- .make_df_sam()

    res <- with_mocked_bindings(
        {
            fRagmentomics:::extract_fragment_features(
                df_sam = df,
                fragment_name = "fragA",
                sample_id = "S_FAIL",
                chr = "chrX", pos = 999L, ref = "G", alt = "C",
                report_bam_info = TRUE,
                report_softclip = TRUE,
                report_5p_3p_bases_fragment = 1L,
                remove_softclip = FALSE,
                fasta_fafile = NULL, fasta_seq = NULL,
                input_mutation_info = "chrX:999:G>C"
            )
        },
        process_fragment_reads_qc = function(df_fragment_reads, chr) "QC_FAIL"
    )

    expect_type(res, "list")
    expect_equal(res$Fragment_QC, "QC_FAIL")
    # a few representative NA placeholders
    expect_true(is.na(res$Fragment_Status_Simple))
    expect_true(is.na(res$Fragment_Size))
    expect_true(is.na(res$Position_5p))
    expect_true(is.na(res$Position_3p))
    expect_true(is.na(res$VAF))
})

test_that("Returns fail object when not exactly one first mate (invalid R1/R2 pairing)", {
    df <- .make_df_sam()

    # Capture the result instead of expecting an error
    result <- with_mocked_bindings(
        fRagmentomics:::extract_fragment_features(
            df_sam = df,
            fragment_name = "fragA",
            sample_id = NA_character_,
            chr = "chr1", pos = 1L, ref = "A", alt = "T",
            report_bam_info = FALSE, report_softclip = FALSE,
            report_5p_3p_bases_fragment = 0L, remove_softclip = FALSE,
            fasta_fafile = NULL, fasta_seq = NULL,
            input_mutation_info = "chr1:1:A>T"
        ),
        # Mock QC to pass the first check
        process_fragment_reads_qc = function(df_fragment_reads, chr) "",
        # Mock flags to simulate invalid R1/R2 pairing
        bamFlagAsBitMatrix = function(flag_vec) {
            cbind(
                isFirstMateRead = c(TRUE, TRUE), # <- invalid (sum != 1)
                isMinusStrand   = c(FALSE, TRUE)
            )
        }
    )

    # Verify that the function returned the correct failure message
    expect_match(result$Fragment_QC, "not a valid R1/R2 pair")

    # Verify that analysis columns are NA (ensuring early return)
    expect_true(is.na(result$Fragment_Status_Simple))
})

test_that("Returns fail object when not exactly one forward and one reverse read", {
    df <- .make_df_sam()

    # Capture the result instead of expecting an error
    result <- with_mocked_bindings(
        fRagmentomics:::extract_fragment_features(
            df_sam = df,
            fragment_name = "fragA",
            sample_id = NA_character_,
            chr = "chr1", pos = 1L, ref = "A", alt = "T",
            report_bam_info = FALSE, report_softclip = FALSE,
            report_5p_3p_bases_fragment = 0L, remove_softclip = FALSE,
            fasta_fafile = NULL, fasta_seq = NULL,
            input_mutation_info = "chr1:1:A>T"
        ),
        # Mock QC to pass the first check
        process_fragment_reads_qc = function(df_fragment_reads, chr) "",
        # Mock flags to simulate invalid strand orientation
        bamFlagAsBitMatrix = function(flag_vec) {
            cbind(
                isFirstMateRead = c(TRUE, FALSE),
                isMinusStrand   = c(FALSE, FALSE) # <- invalid (sum != 1)
            )
        }
    )

    # Verify that the function returned the correct failure message
    expect_match(result$Fragment_QC, "does not have one forward")
    expect_match(result$Fragment_QC, "and one reverse read")

    # Verify that analysis columns are NA
    expect_true(is.na(result$Fragment_Status_Simple))
})

test_that("remove_softclip=TRUE: early return when 5p becomes empty", {
    # use sentinel CIGAR to branch inside our mock
    df <- .make_df_sam(cigar1 = "EMPTY5", cigar2 = "OK")

    res <- with_mocked_bindings(
        {
            fRagmentomics:::extract_fragment_features(
                df_sam = df,
                fragment_name = "fragA",
                sample_id = "Sx",
                chr = "chr1", pos = 123L, ref = "A", alt = "T",
                report_bam_info = TRUE, report_softclip = TRUE,
                report_5p_3p_bases_fragment = 2L, remove_softclip = TRUE,
                fasta_fafile = NULL, fasta_seq = NULL,
                input_mutation_info = "chr1:123:A>T"
            )
        },
        process_fragment_reads_qc = function(df_fragment_reads, chr) "",
        bamFlagAsBitMatrix = function(flag_vec) {
            cbind(isFirstMateRead = c(TRUE, FALSE), isMinusStrand = c(FALSE, TRUE))
        },
        # your trimming function name (whatever it is in your codebase)
        remove_softclip = function(read_stats) {
            if (identical(read_stats$CIGAR, "EMPTY5")) {
                # 5p becomes empty -> triggers early return
                list(SEQ = "", QUAL = "", CIGAR = "", read_length = 0L)
            } else {
                list(
                    SEQ = read_stats$SEQ, QUAL = read_stats$QUAL,
                    CIGAR = read_stats$CIGAR, read_length = nchar(read_stats$SEQ)
                )
            }
        }
    )

    expect_type(res, "list")
    expect_true(res$Fragment_QC != "OK") # early fail after trimming
    expect_true(is.na(res$Fragment_Status_Simple))
    expect_true(is.na(res$Fragment_Size))
})

test_that("remove_softclip=TRUE: early return when 3p becomes empty", {
    df <- .make_df_sam(cigar1 = "OK", cigar2 = "EMPTY3")

    res <- with_mocked_bindings(
        {
            fRagmentomics:::extract_fragment_features(
                df_sam = df,
                fragment_name = "fragA",
                sample_id = "Sy",
                chr = "chr1", pos = 456L, ref = "G", alt = "C",
                report_bam_info = TRUE, report_softclip = TRUE,
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
        }
    )

    expect_type(res, "list")
    expect_true(res$Fragment_QC != "OK")
    expect_true(is.na(res$Fragment_Status_Simple))
    expect_true(is.na(res$Fragment_Size))
})

test_that("TLEN reporting: warning emitted; TLEN is NA_integer_ when missing; Sample_Id present even if NA", {
    df <- .make_df_sam(tlen1 = NA_integer_, tlen2 = NA_integer_)

    expect_warning(
        res <- fRagmentomics:::extract_fragment_features(
            df_sam = df,
            fragment_name = "fragA",
            sample_id = NA_character_,
            chr = "chr1", pos = 888L, ref = "C", alt = "A",
            report_bam_info = TRUE,
            report_softclip = FALSE,
            report_5p_3p_bases_fragment = 0L,
            remove_softclip = FALSE,
            fasta_fafile = NULL, fasta_seq = NULL,
            input_mutation_info = "chr1:888:C>A"
        ),
        regexp = "TLEN is NULL or missing, cannot be reported\\."
    )

    expect_type(res, "list")

    # Sample_Id here -> NA
    expect_true("Sample_Id" %in% names(res))
    expect_true(is.na(res$Sample_Id))

    # TLEN here -> NA_integer_
    expect_true("TLEN" %in% names(res))
    expect_true(is.na(res$TLEN))
    expect_true(is.integer(res$TLEN))

    # No optional output
    expect_false(any(c(
        "Fragment_Bases_5p", "Fragment_Bases_3p",
        "Fragment_Basqs_5p", "Fragment_Basqs_3p"
    ) %in% names(res)))
    expect_false(any(c(
        "Nb_Fragment_Bases_Softclip_5p", "Nb_Fragment_Bases_Softclip_3p"
    ) %in% names(res)))
})


test_that("end_on_reference computes last aligned position using CIGAR operations", {
    # direct unit test for the helper (covers parse_cigar path)
    expect_equal(
        with_mocked_bindings(
            fRagmentomics:::end_on_reference(100L, "7M2D3M"), # ref span = 7 + 2 + 3 = 12
            parse_cigar = function(cigar) {
                data.frame(length = c(7L, 2L, 3L), type = c("M", "D", "M"))
            }
        ),
        111L # 100 + 12 - 1
    )
    # invalid CIGAR or missing POS -> NA
    expect_true(is.na(fRagmentomics:::end_on_reference(NA_integer_, "10M")))
    expect_true(is.na(fRagmentomics:::end_on_reference(100L, "*")))
})

test_that("Check if remove softclip option works properly", {
    df <- .make_df_sam(cigar1 = "5S5M", cigar2 = "5M5S")

    res <- fRagmentomics:::extract_fragment_features(
        df_sam = df,
        fragment_name = "fragA",
        sample_id = "S1",
        chr = "chr1", pos = 12345L, ref = "A", alt = "T",
        report_bam_info = TRUE,
        report_softclip = TRUE,
        report_5p_3p_bases_fragment = 2L,
        remove_softclip = TRUE,
        fasta_fafile = NULL, fasta_seq = NULL,
        input_mutation_info = "chr1:12345:A>T"
    )

    # Function returns a named list (assembled by rbindlist higher up)
    expect_type(res, "list")

    # From original rows: row1 is 5p, row2 is 3p
    expect_equal(res$CIGAR_5p, "5S5M")
    expect_equal(res$CIGAR_3p, "5M5S")

    # TLEN: positive absolute as integer
    expect_equal(res$TLEN, as.integer(abs(.make_df_sam()$TLEN[1])))

    # Optional outputs present
    expect_equal(res$Fragment_Bases_5p, "CG")
    expect_equal(res$Fragment_Bases_3p, "TC")
    expect_equal(res$Nb_Fragment_Bases_Softclip_5p, 0L)
    expect_equal(res$Nb_Fragment_Bases_Softclip_3p, 0L)
})
