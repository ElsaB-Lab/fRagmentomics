test_that("read mut", {
    

    #---------------------------------------
    # chr:pos:ref:alt cases
    #---------------------------------------
    # Case 1 : chr1 
    expect_equal(read_mut("chr1:123:A:T"), data.frame(CHROM = "chr1", POS = 123, REF = "A", ALT = "T", stringsAsFactors = FALSE))
    # Case 2 : Only 1 for chromosome
    expect_equal(read_mut("1:123:A:T"), data.frame(CHROM = "1", POS = 123, REF = "A", ALT = "T", stringsAsFactors = FALSE))
    # Case 3 : X for Chromosome
    expect_equal(read_mut("chrX:456:G:C"), data.frame(CHROM = "chrX", POS = 456, REF = "G", ALT = "C", stringsAsFactors = FALSE))
    # Case 4 : Only X for chromosome
    expect_equal(read_mut("X:789:C:G"), data.frame(CHROM = "X", POS = 789, REF = "C", ALT = "G", stringsAsFactors = FALSE))

    # Special caracters for REF and ALT
    expect_equal(read_mut("1:234:-:A"), data.frame(CHROM = "1", POS = 234, REF = "-", ALT = "A", stringsAsFactors = FALSE))
    expect_equal(read_mut("1:345:.:G"), data.frame(CHROM = "1", POS = 345, REF = ".", ALT = "G", stringsAsFactors = FALSE))
    expect_equal(read_mut("4:111:_:-"), data.frame(CHROM = "4", POS = 111, REF = "_", ALT = "-", stringsAsFactors = FALSE))

    # NA in REF of ALT
    expect_equal(read_mut("2:678:NA:T"), data.frame(CHROM = "2", POS = 678, REF = "NA", ALT = "T", stringsAsFactors = FALSE))
    expect_equal(read_mut("3:910:A:NA"), data.frame(CHROM = "3", POS = 910, REF = "A", ALT = "NA", stringsAsFactors = FALSE))

    # REF or ALT empty
    expect_equal(read_mut("5:222::G"), data.frame(CHROM = "5", POS = 222, REF = "", ALT = "G", stringsAsFactors = FALSE))
    expect_equal(read_mut("6:333:A:"), data.frame(CHROM = "6", POS = 333, REF = "A", ALT = "-", stringsAsFactors = FALSE))
    # Ref and Alt cannot both be empty 
    expect_error(read_mut("7:444::"), 
    "Error: The parameter 'mut' \\(7:444::\\) is not in the expected format \\(.tsv, .vcf, chr:pos:ref:alt\\).")

    # Particular case with "A_" in REF and ALT
    expect_equal(read_mut("8:555:A_:AT"), data.frame(CHROM = "8", POS = 555, REF = "A_", ALT = "AT", stringsAsFactors = FALSE))
    expect_equal(read_mut("8:555:AT:A_"), data.frame(CHROM = "8", POS = 555, REF = "AT", ALT = "A_", stringsAsFactors = FALSE))
    expect_equal(read_mut("8:555:A-:AT"), data.frame(CHROM = "8", POS = 555, REF = "A-", ALT = "AT", stringsAsFactors = FALSE))
    expect_equal(read_mut("8:555:AT:A-"), data.frame(CHROM = "8", POS = 555, REF = "AT", ALT = "A-", stringsAsFactors = FALSE))

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
    read_mut("8:555:TA:TAC,_,,T"),
    data.frame(
        CHROM = rep("8", 4),
        POS = rep(555, 4),
        REF = rep("TA", 4),
        ALT = c("TAC", "_", "", "T"),
        stringsAsFactors = FALSE
        )
    )
    # Invalid case to see eror
    # Miss Ref or Alt
    expect_error(read_mut("chr1:123:A"), 
    "Error: The parameter 'mut' \\(chr1:123:A\\) is not in the expected format \\(.tsv, .vcf, chr:pos:ref:alt\\).")
    # Miss position 
    expect_error(read_mut("chr1:A:T"), 
    "Error: The parameter 'mut' \\(chr1:A:T\\) is not in the expected format \\(.tsv, .vcf, chr:pos:ref:alt\\).")
    # Add extra parameter
    expect_error(read_mut("chr1:123:A:T:extra"), 
    "Error: The parameter 'mut' \\(chr1:123:A:T:extra\\) is not in the expected format \\(.tsv, .vcf, chr:pos:ref:alt\\).")
    # Invalid_format
    expect_error(read_mut("invalid_format"), 
    "Error: The parameter 'mut' \\(invalid_format\\) is not in the expected format \\(.tsv, .vcf, chr:pos:ref:alt\\).")
    # No REF and No ALT
    expect_error(read_mut("1:123::"), 
    "Error: The parameter 'mut' \\(1:123::\\) is not in the expected format \\(.tsv, .vcf, chr:pos:ref:alt\\).")

    #---------------------------------------
    # vcf cases
    #---------------------------------------
    # Case 1: Normal utilisation
    vcf_content <- "##fileformat=VCFv4.2\n##source=Test\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t123\t.\tA\tT\t.\t.\t.\n'2'\t'456'\t.\tG\tC\t454544\tgvgbhj15\t14864_-e_54"
    vcf_file <- tempfile(fileext = ".vcf")
    writeLines(vcf_content, vcf_file)
    
    # Read vcf
    df_vcf <- read_vcf_input(vcf_file)
    
    # Check inside
    expected_df <- data.frame(
        CHROM = c("1", "2"),
        POS = c(123, 456),
        REF = c("A", "G"),
        ALT = c("T", "C"),
        stringsAsFactors = FALSE
    )
    expect_equal(df_vcf, expected_df)
    
    unlink(vcf_file)  # Delete temporary file 

    # Case 2: Does not contain the expected columns
    vcf_content <- "##fileformat=VCFv4.2\n##source=Test\n#CHROM\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t.\tA\tT\t.\t.\t."
    vcf_file <- tempfile(fileext = ".vcf")
    writeLines(vcf_content, vcf_file)

    expect_error(
        read_vcf_input(vcf_file),
        "The VCF file does not contain the expected columns"
    
    )
    
    unlink(vcf_file)

    # Case 3: 0 Row
    vcf_content <- "##fileformat=VCFv4.2\n##source=Test\n#CHROM\tPOS\tREF\tALT"
    vcf_file <- tempfile(fileext = ".vcf")
    writeLines(vcf_content, vcf_file)

    expect_error(
        read_vcf_input(vcf_file),
        "Error: The VCF file does not contain any mutation data."
    )

    unlink(vcf_file)

    # Case 4: Ignore meta data and first line is header
    vcf_content <- "##fileformat=VCFv4.2\n##source=Test\n##comment=ignored\n#CHROM\tPOS\tREF\tALT\nX\t789\tNA\tG\nY\t999\tT\tA"
    vcf_file <- tempfile(fileext = ".vcf")
    writeLines(vcf_content, vcf_file)
    
    df_vcf <- read_vcf_input(vcf_file)
    
    expected_df <- data.frame(
        CHROM = c("X", "Y"),
        POS = c(789, 999),
        REF = c("NA", "T"),
        ALT = c("G", "A"),
        stringsAsFactors = FALSE
    )
    expect_equal(df_vcf, expected_df)
    
    unlink(vcf_file)

    # Case 5: Separator ; instead of tab
    vcf_content <- "##fileformat=VCFv4.2\nCHROM;POS;REF;ALT\n1;123;A;T\n2;456;G;C"
    vcf_file <- tempfile(fileext = ".vcf")
    writeLines(vcf_content, vcf_file)

    expect_error(
        read_vcf_input(vcf_file),
        "Error: The VCF file does not contain the expected columns"
    )

    unlink(vcf_file)

    # Case 6: Multiallelic test
    vcf_content <- "##fileformat=VCFv4.2\n##source=Test\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n8\t555\t.\t-\tT,-,\t.\t.\t."
    vcf_file <- tempfile(fileext = ".vcf")
    writeLines(vcf_content, vcf_file)

    df_vcf <- read_vcf_input(vcf_file)

    expected_df <- data.frame(
        CHROM = rep("8", 3),
        POS = rep(555, 3),
        REF = rep("-", 3),
        ALT = c("T", "-", "-"),
        stringsAsFactors = FALSE
    )

    expect_equal(df_vcf, expected_df)

    unlink(vcf_file)

    # Case 7: Empty lines 
    vcf_content <- "##fileformat=VCFv4.2\n##source=Test\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n8\t555\t.\tA\tT\t.\t.\t.\n\n9\t666\t.\tG\tC\t.\t.\t."
    vcf_file <- tempfile(fileext = ".vcf")
    writeLines(vcf_content, vcf_file)

    df_vcf <- read_vcf_input(vcf_file)

    expect_equal(nrow(df_vcf), 2)

    expected_df <- data.frame(
        CHROM = c("8", "9"),
        POS = c(555, 666),
        REF = c("A", "G"),
        ALT = c("T", "C"),
        stringsAsFactors = FALSE
    )

    expect_equal(df_vcf, expected_df)

    unlink(vcf_file)

    #---------------------------------------
    # tsv cases
    #---------------------------------------
    # Case 1: Correct reading
    tsv_content <- "CHROM\tPOS\tREF\tALT\n1\t123\tA\tT\n2\t456\tG\tC"
    tsv_file <- tempfile(fileext = ".tsv")
    writeLines(tsv_content, tsv_file)
    
    df_tsv <- read_tsv_input(tsv_file)
    
    expected_df <- data.frame(
        CHROM = c("1", "2"),
        POS = c(123, 456),
        REF = c("A", "G"),
        ALT = c("T", "C"),
        stringsAsFactors = FALSE
    )
    expect_equal(df_tsv, expected_df)
    
    unlink(tsv_file)

    # Case 2: Missing column
    tsv_content <- "CHROM\tPOS\tREF\n1\t123\tA\n2\t456\tG"
    tsv_file <- tempfile(fileext = ".tsv")
    writeLines(tsv_content, tsv_file)

    expect_error(
        read_tsv_input(tsv_file),
        "Error: The TSV file does not contain the expected columns"
    )
    
    unlink(tsv_file)

    # Case 3: Empty file
    tsv_content <- "CHROM\tPOS\tREF\tALT" 
    tsv_file <- tempfile(fileext = ".tsv")
    writeLines(tsv_content, tsv_file)

    expect_error(
        read_tsv_input(tsv_file),
        "Error: The TSV file is empty or incorrectly formatted."
    )

    unlink(tsv_file)

    # Case 4: Bad separtor 
    tsv_content <- "CHROM,POS,REF,ALT\n1,123,A,T\n2,456,G,C"
    tsv_file <- tempfile(fileext = ".tsv")
    writeLines(tsv_content, tsv_file)

    expect_error(
        read_tsv_input(tsv_file),
        "Error: The TSV file does not contain the expected columns"
    )

    unlink(tsv_file)

    # Case 5: Empty line ignore 
    tsv_content <- "CHROM\tPOS\tREF\tALT\n1\t123\tA\tT\n\n2\t456\tG\tC"
    tsv_file <- tempfile(fileext = ".tsv")
    writeLines(tsv_content, tsv_file)
    
    df_tsv <- read_tsv_input(tsv_file)

    expect_equal(nrow(df_tsv), 2)

    expected_df <- data.frame(
        CHROM = c("1", "2"),
        POS = c(123, 456),
        REF = c("A", "G"),
        ALT = c("T", "C"),
        stringsAsFactors = FALSE
    )

    expect_equal(df_tsv, expected_df)

    unlink(tsv_file)

    # Case 6: Mutliallelic test
    tsv_content <- "CHROM\tPOS\tREF\tALT\n8\t555\tA\tT,-,\n9\t666\tG\tC,A"
    tsv_file <- tempfile(fileext = ".tsv")
    writeLines(tsv_content, tsv_file)

    df_tsv <- read_tsv_input(tsv_file)

    expected_df <- data.frame(
        CHROM = c("8", "8", "8", "9", "9"),
        POS = c(555, 555, 555, 666, 666),
        REF = c("A", "A", "A", "G", "G"),
        ALT = c("T", "-", "-", "C", "A"),
        stringsAsFactors = FALSE
    )

    expect_equal(df_tsv, expected_df)

    unlink(tsv_file)

    # Case 7: Logical vector with tsv
    tsv_content <- "CHROM\tPOS\tREF\tALT\n8\t555\tA\tT\n9\t666\tG\tNA"
    tsv_file <- tempfile(fileext = ".tsv")
    writeLines(tsv_content, tsv_file)

    df_tsv <- read_tsv_input(tsv_file)

    expected_df <- data.frame(
        CHROM = c("8", "9"),
        POS = c(555, 666),
        REF = c("A", "G"),
        ALT = c("T", "NA"),
        stringsAsFactors = FALSE
    )

    expect_equal(df_tsv, expected_df)

    unlink(tsv_file)
    }
)