# Project : ElsaBLab_fRagmentomics

#' Parser VCF
#'
#' This function checks if the VCF file contains the expected structure and columns "CHROM", "POS", "REF", "ALT".
#'
#' @param vcf_file Path to the VCF file.
#' @return A cleaned data frame with required columns.
#' 
#' @noRd
parser_vcf <- function(vcf_file) {
    # Determine if the file is compressed
    is_compressed <- grepl("\\.vcf\\.gz$", vcf_file)

    # Read all lines while handling compressed files
    con <- if (is_compressed) gzfile(vcf_file, "rt") else vcf_file
    lines <- readLines(con)
    close(con)

    # Keep only the useful info (header + data)
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

    # Convert POS column into int (Force to NA if str)
    vcf_subset[["POS"]] <- suppressWarnings(as.integer(vcf_subset[["POS"]]))

    # Check if there is NA in POS column after conversion
    invalid_pos_rows <- vcf_subset[is.na(vcf_subset[["POS"]]), ]

    # If line have no convertible POS value into int, print it
    if (nrow(invalid_pos_rows) > 0) {
        warning("Warning: The following rows have non-integer values in the POS column:")
        print(invalid_pos_rows)

    # Remove empty lines
    vcf_subset <- vcf_subset[rowSums(is.na(vcf_subset) | vcf_subset == "") != ncol(vcf_subset), ]

    if (nrow(vcf_subset) == 0) {
        stop("Error: The VCF file does not contain any mutation data.")
    }

    return(vcf_subset)
    }
}


#' Parser TSV
#'
#' This function checks if the TSV file contains the expected structure and columns.
#'
#' @param tsv_file Path to the TSV file.
#' @return A cleaned data frame with required columns.
#' 
#' @noRd
parser_tsv <- function(tsv_file) {
    # Determine if the file is compressed
    is_compressed <- grepl("\\.tsv\\.gz$", tsv_file)

    # Read the file while handling compression
    con <- if (is_compressed) gzfile(tsv_file, "rt") else tsv_file
    df_tsv <- tryCatch(
        read.table(
            con,
            header = TRUE,
            sep = "\t",
            colClasses = "character",    # Force everything to character
            na.strings = c(),            # Prevent automatic conversion of "NA" to <NA>
            stringsAsFactors = FALSE
        ),
        error = function(e) {
            stop("Error: Failed to read the TSV file. Ensure it's properly formatted.")
        }
    )
    if (is_compressed) close(con)

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

    # Convert POS column into int (Force to NA if str)
    tsv_subset[["POS"]] <- suppressWarnings(as.integer(tsv_subset[["POS"]]))

    # Check if there is NA in POS column after conversion
    invalid_pos_rows <- tsv_subset[is.na(tsv_subset[["POS"]]), ]

    # If line have no convertible POS value into int, print it
    if (nrow(invalid_pos_rows) > 0) {
        warning("Warning: The following rows have non-integer values in the POS column:")
        print(invalid_pos_rows)

    # Remove empty rows
    tsv_subset <- tsv_subset[rowSums(is.na(tsv_subset) | tsv_subset == "") != ncol(tsv_subset), ]

    # Check if the TSV has data
    if (nrow(tsv_subset) < 2) {
        stop("Error: The TSV file does not contain any mutation data.")
    }

    return(tsv_subset)
    }
}

#' Parser chr:pos:ref:alt
#'
#' This function checks if the mutation info contains the expected well structured informations.
#'
#' @param mut String chr:pos:ref:alt
#' @return A cleaned data frame with required informations
#' 
#' @noRd
parser_mut_str <- function(mut) {
    # If mut is finished by ":", the programm will add - at the end to use strsplit
    if (grepl(":$", mut)) {
      mut <- paste0(mut, "-")
    }

    # Devide the 4 parts of the mutation expected
    parts <- unlist(strsplit(mut, ":"))

    # Check if mutation contains exactly these parts (CHR, POS, REF, ALT)
    if (length(parts) != 4) {
      stop(paste0("Error: The parameter 'mut' (", mut, ") is not in the expected format (.tsv, .vcf, chr:pos:ref:alt)."))
    }

    chr <- parts[1]

    # Only take into account integer
    if (!grepl("^[0-9]+$", pos_str)) {
        pos <- NA
    } else {
        pos <- suppressWarnings(as.integer(parts[2]))
    }

    ref <- ref <- parts[3]
    alt <- alt <- parts[4]

    # Return the info into a list 
    mut_df <- data.frame(CHROM = chr, POS = pos, REF = ref, ALT = alt, stringsAsFactors = FALSE)
}