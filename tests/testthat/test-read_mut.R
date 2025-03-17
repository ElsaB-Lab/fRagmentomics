test_that("read mut", {
    
    vcf_test <- load(system.file("testdata", "mutations_test.vcf", package="fRagmentomics"))
    tsv_test <- load(system.file("testdata", "mutations_test.tsv", package="fRagmentomics"))
    vcf_test_compressed <- load(system.file("testdata", "mutations_test.vcf.gz", package="fRagmentomics"))
    tsv_test_compressed <- load(system.file("testdata", "mutations_test.tsv.gz", package="fRagmentomics"))


    #---------------------------------------
    # chr:pos:ref:alt cases
    #---------------------------------------
    # Normal cases
    expect_equal(read_mut("chr1:123:A:T"), data.frame(CHROM = "chr1", POS = 123, REF = "A", ALT = "T", stringsAsFactors = FALSE))
    # X for Chromosome
    expect_equal(read_mut("chrX:456:G:C"), data.frame(CHROM = "chrX", POS = 456, REF = "G", ALT = "C", stringsAsFactors = FALSE))

    # Special caracters for REF and ALT
    expect_equal(read_mut("4:111:_:-"), data.frame(CHROM = "4", POS = 111, REF = "_", ALT = "-", stringsAsFactors = FALSE))

    # NA in REF of ALT
    expect_equal(read_mut("2:678:NA:T"), data.frame(CHROM = "2", POS = 678, REF = "NA", ALT = "T", stringsAsFactors = FALSE))

    # Mutliallelic cases 
    expect_equal(
    read_mut("8:555:A:T,-,"),
    data.frame(
        CHROM = rep("8", 3),
        POS = rep(555, 3),
        REF = rep("A", 3),
        ALT = c("T", "-", "-"),
        stringsAsFactors = FALSE
        )
    )
    expect_equal(
    read_mut("8:555:TAC,_,,T:TA"),
    data.frame(
        CHROM = rep("8", 4),
        POS = rep(555, 4),
        REF = c("TAC", "_", "", "T"),
        ALT = rep("TA", 4),
        stringsAsFactors = FALSE
        )
    )
    # Invalid case to see eror
    # Miss Ref or Alt
    expect_error(read_mut("chr1:123:A"), 
    "Error: The parameter 'mut' \\(chr1:123:A\\) is not in the expected format \\(.tsv, .vcf, chr:pos:ref:alt\\).")
    # Add extra parameter
    expect_error(read_mut("chr1:123:A:T:extra"), 
    "Error: The parameter 'mut' \\(chr1:123:A:T:extra\\) is not in the expected format \\(.tsv, .vcf, chr:pos:ref:alt\\).")
    # Invalid_format
    expect_error(read_mut("invalid_format"), 
    "Error: The parameter 'mut' \\(invalid_format\\) is not in the expected format \\(.tsv, .vcf, chr:pos:ref:alt\\).")

    #---------------------------------------
    # vcf cases
    #---------------------------------------
    # Case 1: VCF 
    df_vcf_test <- read_mut(vcf_test)
    print("df_vcf_test:")
    print(df_vcf_test)

    expected_results <- data.frame(
        CHROM = c("chr1", "2", "chr3", "chr22", "22", "chrX", "chrY", "chr4", "chr5", "chr6",
                    "chr7", "chr7", "chr8", "chr8", "chrX", "chrY", "42453", "."),
        POS = c(12345, 67890, 101112, 56789, 56789, 54321, 99999, 1234, 5678, 91011,
                121314, 121314, 151617, 151617, "NA", "NA", 181920, 181920),
        REF = c("A", "GT", "A.", "NA", ".", "A", ".", "-", ".", "G",
                "T", "T", "G", "G", "A", "AC", "AC"),
        ALT = c("T", "C", "A-", "NA", "-", ".", "C", "A", ".", "NA",
                "C", "C", "A", "T", "-", "C", "-", "-"),
        stringsAsFactors = FALSE
    )

    expect_equal(df_vcf_test, expected_results)

    # Case 2: Compressed VCF
    df_vcf_test2 <- read_mut(vcf_test_compressed)
    print("df_vcf_test2:")
    print(df_vcf_test2)

    expected_results2 <- data.frame(
        CHROM = c("chr1", "2", "chr3", "chr22", "22", "chrX", "chrY", "chr4", "chr5", "chr6",
                    "chr7", "chr7", "chr8", "chr8", "chrX", "chrY", "42453", "."),
        POS = c(12345, 67890, 101112, 56789, 56789, 54321, 99999, 1234, 5678, 91011,
                121314, 121314, 151617, 151617, "NA", "NA", 181920, 181920),
        REF = c("A", "GT", "A.", "NA", ".", "A", ".", "-", ".", "G",
                "T", "T", "G", "G", "A", "AC", "AC"),
        ALT = c("T", "C", "A-", "NA", "-", ".", "C", "A", ".", "NA",
                "C", "C", "A", "T", "-", "C", "-", "-"),
        stringsAsFactors = FALSE
    )

    expect_equal(df_vcf_test2, expected_results2)

    #---------------------------------------
    # tsv cases
    #---------------------------------------
    # Case 1: TSV
    df_tsv_test <- read_mut(tsv_test)
    print("df_tsv_test:")
    print(df_tsv_test)

    expected_results3 <- data.frame(
        CHROM = c("chr1", "2", "chr3", "chr22", "22", "chrX", "chrY", "chr4", "chr5", "chr6",
                    "chr7", "chr7", "chr8", "chr8", "chrX", "chrY", "42453", "."),
        POS = c(12345, 67890, 101112, 56789, 56789, 54321, 99999, 1234, 5678, 91011,
                121314, 121314, 151617, 151617, "NA", "NA", 181920, 181920),
        REF = c("A", "GT", "A.", "NA", ".", "A", ".", "-", ".", "G",
                "T", "T", "G", "G", "A", "AC", "AC"),
        ALT = c("T", "C", "A-", "NA", "-", ".", "C", "A", ".", "NA",
                "C", "C", "A", "T", "-", "C", "-", "-"),
        stringsAsFactors = FALSE
    )

    expect_equal(df_tsv_test, expected_results3)

    # Case 2: Compressed TSV
    df_tsv_test2 <- read_mut(tsv_test_compressed)
    print("df_tsv_test2:")
    print(df_tsv_test2)

    expected_results4 <- data.frame(
        CHROM = c("chr1", "2", "chr3", "chr22", "22", "chrX", "chrY", "chr4", "chr5", "chr6",
                    "chr7", "chr7", "chr8", "chr8", "chrX", "chrY", "42453", "-"),
        POS = c(12345, 67890, 101112, 56789, 56789, 54321, 99999, 1234, 5678, 91011,
                121314, 121314, 151617, 151617, "NA", "NA", 181920, 181920),
        REF = c("A", "GT", "A.", "NA", ".", "A", ".", "-", ".", "G",
                "T", "T", "G", "G", "A", "AC", "AC"),
        ALT = c("T", "C", "A-", "NA", "-", ".", "C", "A", ".", "NA",
                "C", "C", "A", "T", "-", "C", "-", "-"),
        stringsAsFactors = FALSE
    )

    expect_equal(df_tsv_test2, expected_results4)
    }
)