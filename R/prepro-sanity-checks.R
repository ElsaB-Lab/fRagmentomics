# Project : ElsaBLab_fRagmentomics


#' Sanity Check VCF
#'
#' This function checks if the VCF file contains the expected structure and columns "CHROM", "POS", "REF", "ALT".
#'
#' @param vcf_file Path to the VCF file.
#' @return A cleaned data frame with required columns.
#' 
#' @noRd
sanity_check_vcf <- function(vcf_file) {
    # Check if file exists
    if (!file.exists(vcf_file)) {
        stop(paste("Error: The file", vcf_file, "does not exist."))
    }

    # Read all files and keep only the usefull info (header + data)
    lines <- readLines(vcf_file)
    data_lines <- lines[!grepl("^##", lines)]

    # To do the read.table, we need at list one row of data
    if (length(data_lines) < 2) {
        stop("Error: The VCF file does not contain any mutation data.")
    }

    # Get the header to create the final file and the expected columns
    header_line <- data_lines[1]
    header <- strsplit(header_line, "	")[[1]]
    header[1] <- sub("^#", "", header[1])
    expected_cols <- c("CHROM", "POS", "REF", "ALT")

    # Read the data while avoiding auto-conversion of 'T' or 'NA' to logical/NA
    df_vcf <- read.table(
        vcf_file,
        header = FALSE,
        comment.char = "#",
        sep = "\t",
        colClasses = "character",    # Force everything to character
        na.strings = c(),            # Prevent automatic conversion of "NA" to <NA>
        stringsAsFactors = FALSE
    )

    # Check if the read.table worked propely
    if (length(header) != ncol(df_vcf)) {
        stop("Error: The number of columns in the VCF file does not match the header line.")
    }

    # Define the header
    colnames(df_vcf) <- header

    # Check the columns of the final df
    if (!all(expected_cols %in% colnames(df_vcf))) {
        stop(paste("Error: The VCF file does not contain the expected columns (CHROM, POS, REF, ALT).",
                   "Found columns:", paste(colnames(df_vcf), collapse = ", ")))
    }

    # Create the df with only the necessary column
    vcf_subset <- df_vcf[, expected_cols, drop = FALSE]

    # Convert POS to integer explicitly
    # If there's any invalid integer, we'll catch it via `is.na(...)`
    vcf_subset[["POS"]] <- as.integer(vcf_subset[["POS"]])
    if (any(is.na(vcf_subset[["POS"]]))) {
        stop("Error: The POS column contains non-integer values.")
    }

    # Remove empty lines
    vcf_subset <- vcf_subset[rowSums(is.na(vcf_subset) | vcf_subset == "") != ncol(vcf_subset), ]

    if (nrow(vcf_subset) == 0) {
        stop("Error: The VCF file does not contain any mutation data.")
    }

    return(vcf_subset)
}


#' Sanity Check TSV
#'
#' This function checks if the TSV file contains the expected structure and columns.
#'
#' @param tsv_file Path to the TSV file.
#' @return A cleaned data frame with required columns.
sanity_check_tsv <- function(tsv_file) {
    # Check if file exists
    if (!file.exists(tsv_file)) {
        stop(paste("Error: The file", tsv_file, "does not exist."))
    }

    # Read the data while avoiding auto-conversion of 'T' or 'NA' to logical/NA
    df_tsv <- read.table(
        tsv_file,
        header = TRUE,
        sep = "\t",
        colClasses = "character",    # Force everything to character
        na.strings = c(),            # Prevent automatic conversion of "NA" to <NA>
        stringsAsFactors = FALSE
    )

    # Check if the file is empty
    if (nrow(df_tsv) == 0) {
        stop("Error: The TSV file is empty or incorrectly formatted.")
    }

    # Define expected columns
    expected_cols <- c("CHROM", "POS", "REF", "ALT")

    # Ensure all expected columns exist
    if (!all(expected_cols %in% colnames(df_tsv))) {
        stop(paste("Error: The TSV file does not contain the expected columns (CHROM, POS, REF, ALT).",
                   "Found columns:", paste(colnames(df_tsv), collapse = ", ")))
    }

    # Extract only the necessary columns
    tsv_subset <- df_tsv[, expected_cols, drop = FALSE]

    # Convert POS to integer explicitly
    # If there's any invalid integer, we'll catch it via `is.na(...)`
    tsv_subset[["POS"]] <- as.integer(tsv_subset[["POS"]])
    if (any(is.na(tsv_subset[["POS"]]))) {
        stop("Error: The POS column contains non-integer values.")
    }

    # Remove empty rows
    tsv_subset <- tsv_subset[rowSums(is.na(tsv_subset) | tsv_subset == "") != ncol(tsv_subset), ]

    # Check if the TSV has data
    if (nrow(tsv_subset) < 2) {
        stop("Error: The TSV file does not contain any mutation data.")
    }

    return(tsv_subset)
}