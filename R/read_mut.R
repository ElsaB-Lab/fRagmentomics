#' Read variant information from multiple sources
#'
#' @description This function acts as a dispatcher to read variant information from a file (VCF or TSV) or a character
#' string. It automatically detects the input format and uses the appropriate parser. All multi-allelic sites are split
#' into separate, bi-allelic rows.
#'
#' @inheritParams run_fRagmentomics
#'
#' @return A 'data.frame' where each row represents a single bi-allelic variant, with the columns 'CHROM', 'POS', 'REF', and 'ALT'.
#'
#' @keywords internal
read_mut <- function(mut) {
  if (grepl("\\.tsv(.gz)?$", mut)) {
    df_mut <- read_tsv_input(mut)
  } else if (grepl("\\.vcf(.gz)?$", mut)) {
    df_mut <- read_vcf_input(mut)
  } else if (grepl("^[^:]+:[^:]+:[^:]+:[^:]+$", mut)) {
    df_mut <- read_str_input(mut)
  } else {
    stop(sprintf(
      "The parameter 'mut' ('%s') is not in the expected format (.tsv, .vcf, chr:pos:ref:alt).",
      mut
    ))
  }

  df_mut
}

#' Read and process variants from a VCF file
#'
#' @description This function reads variant data from a VCF file, uses a low-level parser for initial validation
#' and then expands any multi-allelic records into a standard bi-allelic format.
#'
#' @param vcf_file Path to the VCF file.
#'
#' @return A processed dataframe with normalized variants.
#'
#' @noRd
read_vcf_input <- function(vcf_file) {
  vcf_subset <- parser_vcf(vcf_file)
  # Create a row for each alt in multiallelics cases
  vcf_subset_without_multiallel <- expand_multiallelics(vcf_subset)
  vcf_subset_without_multiallel
}

#' Read and process variants from a TSV file
#'
#' @description This function reads variant data from a TSV file, uses a low-level parser for initial validation and
#' then expands any multi-allelic records into a standard bi-allelic format.
#'
#' @param tsv_file Path to the TSV file.
#'
#' @return A processed dataframe with normalized variants.
#'
#' @noRd
read_tsv_input <- function(tsv_file) {
  tsv_subset <- parser_tsv(tsv_file)
  tsv_subset_without_multiallel <- expand_multiallelics(tsv_subset)
  tsv_subset_without_multiallel
}


#' Read and process a variant from a string
#'
#' @description This function reads variant data from a formatted character string, uses a low-level parser for validation
#' and then expands any multi-allelic records into a standard bi-allelic format.
#'
#' @param mut Infos chr:pos:ref:alt
#'
#' @return A processed dataframe with normalized variants.
#'
#' @noRd
read_str_input <- function(mut) {
  mut_df <- parser_str(mut)
  mut_df_without_multiallel <- expand_multiallelics(mut_df)
  mut_df_without_multiallel
}


#' Parse a VCF file for variant information
#'
#' @description This low-level parser reads a VCF file, performs basic
#' structural checks, and extracts the essential 'CHROM', 'POS', 'REF', and
#' 'ALT' columns. It does not expand multi-allelic sites.
#'
#' @param vcf_file Path to the VCF file.
#'
#' @return A cleaned dataframe with required columns.
#'
#' @importFrom utils read.table
#'
#' @noRd
parser_vcf <- function(vcf_file) {
  # Determine if the file is compressed
  is_compressed <- grepl("\\.vcf\\.gz$", vcf_file)

  # Read all lines while handling compressed files
  con <- if (is_compressed) gzfile(vcf_file, "rt") else vcf_file
  lines <- tryCatch(
    readLines(con),
    warning = function(w) {
      # If the warning is the specific one about an incomplete final line,
      # ignore it by invoking a restart that muffles the warning.
      if (grepl("incomplete final line found", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      } else {
        # For any other warning, re-issue it as it might be important.
        warning(w)
      }
    }
  )
  if (inherits(con, "connection")) close(con)

  # Keep only the lines that are not metadata (##)
  data_lines <- lines[!grepl("^##", lines)]

  # A valid VCF must have at least a header (#) and one data row
  if (length(data_lines) < 2) {
    stop("The VCF file does not contain any mutation data.")
  }

  # Extract the header from the first non-metadata line
  header_line <- data_lines[1]
  header <- strsplit(header_line, "\t")[[1]]
  header[1] <- sub("^#", "", header[1])

  # Read the data while avoiding auto-conversion of 'T' or 'NA' to logical/NA
  df_vcf <- read.table(
    vcf_file,
    header = FALSE,
    comment.char = "#",
    sep = "\t",
    col.names = header,
    colClasses = "character", # Force everything to character
    na.strings = "", # Prevent automatic conversion of "NA" to <NA>
    stringsAsFactors = FALSE,
    fill = TRUE
  )

  # Keep only the essential columns
  df_vcf <- df_vcf[, c("CHROM", "POS", "REF", "ALT"), drop = FALSE]

  # Check if the expected columns were found
  if (ncol(df_vcf) != 4) {
    stop(sprintf(
      "The VCF file does not seem to have the required columns (CHROM, POS, REF, ALT). Found columns: %s",
      paste(colnames(df_vcf), collapse = ", ")
    ))
  }

  # Identify rows with non-numeric values in the POS column before conversion
  invalid_pos_mask <- !grepl("^[0-9]+$", df_vcf[["POS"]]) & !is.na(df_vcf[["POS"]])

  # If any are found, issue a single, clear warning
  if (any(invalid_pos_mask)) {
    invalid_rows <- which(invalid_pos_mask)
    warning(sprintf(
      "The following rows have non-integer values in the POS column and will be converted to NA: %s",
      paste(invalid_rows, collapse = ", ")
    ))
  }

  # Convert the POS column to integer
  df_vcf[["POS"]] <- as.integer(df_vcf[["POS"]])

  # Remove any rows that are entirely empty or NA
  df_vcf <- df_vcf[rowSums(is.na(df_vcf) | df_vcf == "") != ncol(df_vcf), ]

  if (nrow(df_vcf) == 0) {
    stop("The VCF file does not contain any mutation data after parsing.")
  }

  df_vcf
}


#' Parse a TSV file for variant information
#'
#' @description This low-level parser reads a TSV file, performs basic structural checks, and extracts the essential
#' 'CHROM', 'POS', 'REF', and 'ALT' columns. It does not expand multi-allelic sites.
#'
#' @param tsv_file Path to the TSV file.
#'
#' @return A cleaned dataframe with required columns.
#'
#' @importFrom readr read_tsv
#'
#' @noRd
parser_tsv <- function(tsv_file) {
  # Read the file while handling compression
  df_tsv <- tryCatch(
    {
      readr::read_tsv(
        file = tsv_file,
        col_types = readr::cols_only(
          CHROM = readr::col_character(),
          POS   = readr::col_character(),
          REF   = readr::col_character(),
          ALT   = readr::col_character()
        ),
        na = character(),
        progress = FALSE
      )
    },
    error = function(e) {
      stop("Failed to read the TSV file. Ensure it's properly formatted.")
    }
  )

  # Check if the read.table worked propely
  if (ncol(df_tsv) != 4) {
    stop(sprintf(
      "The number of columns in the TSV file does not match the expected columns CHR POS REF ALT. Found columns: %s",
      paste(colnames(df_tsv), collapse = ", ")
    ))
  }

  # Identify lines with non numeric values in POS
  invalid_pos_mask <- !grepl("^[0-9]+$", df_tsv[["POS"]]) & !is.na(df_tsv[["POS"]])

  # Warning
  if (any(invalid_pos_mask)) {
    invalid_rows <- which(invalid_pos_mask)
    warning(sprintf(
      "The following rows have non-integer values in the POS column and will be converted to NA: %s",
      paste(invalid_rows, collapse = ", ")
    ))
  }

  # Conversion of the position in integer
  df_tsv[["POS"]] <- as.integer(df_tsv[["POS"]])

  # Remove empty rows
  df_tsv <- df_tsv[rowSums(is.na(df_tsv) | df_tsv == "") != ncol(df_tsv), ]

  # Check if the TSV has data
  if (nrow(df_tsv) == 0) {
    stop("The TSV file does not contain any mutation data.")
  }

  df_tsv
}


#' Parse a string for variant information
#'
#' @description This low-level parser reads a "chr:pos:ref:alt" formatted string, performs basic structural checks,
#' and extracts the four components into a data frame. It does not expand multi-allelic sites.
#'
#' @param mut String chr:pos:ref:alt
#'
#' @return A cleaned dataframe with required informations
#'
#' @noRd
parser_str <- function(mut) {
  # If mut is finished by ":", the programm will add - at the end to use strsplit
  if (grepl(":$", mut)) {
    mut <- paste0(mut, "-")
  }

  # Devide the 4 parts of the mutation expected
  parts <- unlist(strsplit(mut, ":"))

  # Check if mutation contains exactly these parts (CHR, POS, REF, ALT)
  if (length(parts) != 4) {
    stop(sprintf(
      "The parameter 'mut' ('%s') is not in the expected format (.tsv, .vcf, chr:pos:ref:alt).",
      mut
    ))
  }

  chr <- parts[1]

  # Only take into account integer
  if (!grepl("^[0-9]+$", parts[2])) {
    warning(sprintf("Position value '%s' is not a valid integer and will be treated as NA.", parts[2]))
    pos <- NA
  } else {
    pos <- as.integer(parts[2])
  }

  ref <- ref <- parts[3]
  alt <- alt <- parts[4]

  # Return the info into a list
  mut_df <- data.frame(CHROM = chr, POS = pos, REF = ref, ALT = alt, stringsAsFactors = FALSE)

  # Check if the TSV has data
  if (nrow(mut_df) == 0) {
    stop("The mutation string does not contain any mutation data.")
  }

  mut_df
}


#' Expand multi-allelic variant records into bi-allelic format
#'
#' @description This function takes a data frame of variants and expands rows where the 'ALT' column contains multiple
#' comma-separated alleles. A new row is created for each alternate allele, while the 'CHROM', 'POS', and 'REF' values are
#' duplicated.
#'
#' @param df A dataframe with columns CHROM, POS, REF, and ALT.
#'
#' @return A dataframe where each row contains only a single alternative allele.
#'
#' @importFrom tidyr separate_rows
#' @noRd
expand_multiallelics <- function(df) {
  # Identify rows where REF contains a comma
  is_multi_ref <- grepl(",", df$REF)

  # If any such rows exist, issue a warning and remove them
  if (any(is_multi_ref)) {
    warning("REF can not be multiallelic. These variants have been removed.")
    df <- df[!is_multi_ref, ]
  }

  # If no data remains after filtering, stop
  if (nrow(df) == 0) {
    stop("After reading and filtering, the mutation information is empty. No valid mutation data found.")
  }

  # Replace NA or "" with "-"
  df$ALT[is.na(df$ALT) | df$ALT == ""] <- "-"
  # Append a "-" if the string ends with a comma (e.g., "A," becomes "A,-")
  df$ALT <- sub(",$", ",-", df$ALT)

  # Split alternative alleles and expand the dataframe
  expanded_df <- tidyr::separate_rows(df, ALT, sep = ",")

  # Return the df
  return(as.data.frame(expanded_df))
}
