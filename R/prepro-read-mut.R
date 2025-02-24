read_vcf <- function (vcf_file) {
    # Open the file and skip the metadata
    con <- file(vcf_file, "r")
    repeat {
        header_line <- readLines(con, n = 1)
        # We stop at the first line which doesn't start by ##
        if (!grepl("^##", header_line)) {
            break
        }
    }
    close(con)

    # Define the header 
    header <- strsplit(header_line, "\t")[[1]]
    header[1] <- sub("^#", "", header[1])

    # Define expected columns 
    expected_cols <- c("CHROM", "POS", "REF", "ALT")

    # Check if the VCF contains the good columns
    if (!all(expected_cols %in% header)) {
        stop(paste("Error: The VCF file does not contain the expected columns (CHROM, POS, REF, ALT).",
                   "Found columns:", paste(header, collapse = ", ")))
    }

    # Find the indices of the columns to extract
    col_indices <- match(expected_cols, header)

    # Read the VCF file 
    # Delete all the lines started by # exept the one just before the data
    df_vcf <- read.table(vcf_file, header = FALSE, comment.char = "#", sep = "\t", stringsAsFactors = FALSE)

    # Check if the vcf has rows 
    if (nrow(df_vcf) == 0) {
    stop("Error: The VCF file does not contain any mutation data.")
    }
    
    # We only keep the columns chr, pos, ref, alt from the vcf
    vcf_subset <- df_vcf[, col_indices, drop = FALSE]
    return(vcf_subset)
}

read_tsv <- function(tsv_file) {
    # Read the TSV file with header
    df_tsv <- read.table(tsv_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # Define expected columns
    expected_cols <- c("CHROM", "POS", "REF", "ALT")

    # Check if the TSV contains the expected columns
    if (!all(expected_cols %in% colnames(df_tsv))) {
        stop(paste("Error: The TSV file does not contain the expected columns (CHROM, POS, REF, ALT).",
                   "Found columns:", paste(colnames(df_tsv), collapse = ", ")))
    }

    # Find the indices of the required columns
    col_indices <- match(expected_cols, colnames(df_tsv))

    # Extract only the necessary columns
    tsv_subset <- df_tsv[, col_indices, drop = FALSE]

    # Check if the TSV has data
    if (nrow(tsv_subset) == 0) {
        stop("Error: The TSV file does not contain any mutation data.")
    }

    return(tsv_subset)
}

read_mut <- function(mut) {
  if (grepl("\\.tsv$", mut)) {

    # Return a datafram without header, with all the mutations from the vcf
    # Columns : chr pos ref alt
    df_mut <- read_tcv(mut)
    return(df_mut)

  } else if (grepl("\\.vcf$", mut)) {
    
    # Return a datafram without header, with all the mutations from the vcf
    # Columns : chr pos ref alt
    df_mut <- read_vcf(mut)
    return(df_mut)

  } else if (grepl("^chr[0-9XY]+:[0-9]+:[ACGT]+:[ACGT]+$", mut)) {

    # Extract chr pos ref alt
    parts <- unlist(strsplit(mut, ":"))

    chr <- parts[1]
    pos <- pos <- as.integer(parts[2])
    ref <- ref <- parts[3]
    alt <- alt <- parts[4]

    # Retunr the info into a list 
    return(data.frame(CHROM = chr, POS = pos, REF = ref, ALT = alt, stringsAsFactors = FALSE))

  } else {
    stop("Error : The parameter 'mut' is not in the expected format (.tsv, .vcf, chr:pos:ref:alt).")
    }
}


