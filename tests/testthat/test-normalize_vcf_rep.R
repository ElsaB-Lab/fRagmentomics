test_that("normalize_vcf_rep", {
    fasta_38 <- "/mnt/beegfs02/database/bioinfo/Index_DB/BWA/0.7.17/UCSC/hg38/hg38.fa"

    # Valid cases
    expect_equal(
        normalize_vcf_rep("chr17", 7676126, "TG", "TGT", fasta_38),
        data.frame(chr = "chr17", pos = 7676127, ref = "G", alt = "GT", stringsAsFactors = FALSE)
    )

    expect_equal(
        normalize_vcf_rep("chr17", 7676126, "TGTA", "TG", fasta_38),
        data.frame(chr = "chr17", pos = 7676127, ref = "GTA", alt = "G", stringsAsFactors = FALSE)
    )

    expect_equal(
        normalize_vcf_rep("chr17", 7676126, "TGT", "TGTAGG", fasta_38),
        data.frame(chr = "chr17", pos = 7676128, ref = "T", alt = "TAGG", stringsAsFactors = FALSE)
    )

    expect_equal(
        normalize_vcf_rep("chr17", 7676126, "TGTA", "AG", fasta_38),
        data.frame(
        chr = c("chr17", "chr17"), 
        pos = c(7676126, 7676127), 
        ref = c("T", "GTA"), 
        alt = c("A", "G"),
        stringsAsFactors = FALSE
        )
    )

    expect_equal(
        normalize_vcf_rep("chr17", 7676125, "GTGT", "GTATCC", fasta_38),
        data.frame(
        chr = c("chr17", "chr17"), 
        pos = c(7676127, 7676128), 
        ref = c("G", "T"), 
        alt = c("A", "TCC"),
        stringsAsFactors = FALSE
        )
    )

    expect_equal(
        normalize_vcf_rep("chr17", 7676126, "TG", "TG", fasta_38),
        data.frame(chr = "chr17", pos = 7676126, ref = "TG", alt = "TG", stringsAsFactors = FALSE)
    )

    # Error cases
    expect_error(normalize_vcf_rep("chr17", 7676126, "AGT", "A", fasta_38))
    expect_error(normalize_vcf_rep("chr17", 7676126, "TAT", "TA", fasta_38))
    expect_error(normalize_vcf_rep("chr17", 7676126, "T", "", fasta_38))
    expect_error(normalize_vcf_rep("chr17", 7676126, "T", ".", fasta_38))
    expect_error(normalize_vcf_rep("chr17", 7676126, "T", "-", fasta_38))
    expect_error(normalize_vcf_rep("chr17", 7676126, "", "GG", fasta_38))
    expect_error(normalize_vcf_rep("chr17", 7676126, ".", "GG", fasta_38))
    expect_error(normalize_vcf_rep("chr17", 7676126, "-", "GG", fasta_38))

    }
)