# Project : ElsaBLab_fRagmentomics

#' Parser VCF
#'
#' @param vcf_file Path to the VCF file.
#'
#' @return A cleaned data frame with required columns.
#'
#' @importFrom utils read.table
#'
#' @noRd
parser_vcf <- function(vcf_file) {
  # Determine if the file is compressed
  is_compressed <- grepl("\\.vcf\\.gz$", vcf_file)

  # Read all lines while handling compressed files
  con <- if (is_compressed) gzfile(vcf_file, "rt") else vcf_file
  lines <- suppressWarnings(readLines(con))
  if (inherits(con, "connection")) close(con)

  # Keep only the useful info (header + data)
  data_lines <- lines[!grepl("^##", lines)]

  # To do the read.table, we need at list one row of data
  if (length(data_lines) < 2) {
    stop("Error: The VCF file does not contain any mutation data.")
  }

  # Get the header to create the final file and the expected columns
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

  # Name the columns
  df_vcf <- df_vcf[, c("CHROM", "POS", "REF", "ALT"), drop = FALSE]

  # Check if the read.table worked propely
  if (ncol(df_vcf) != 4) {
    stop(paste(
      "Error: The number of columns in the VCF file does not match the expected columns CHR POS REF ALT.",
      "Found columns:", paste(colnames(df_vcf), collapse = ", ")
    ))
  }

  # Convert POS column into int (Force to NA if str)
  df_vcf[["POS"]] <- suppressWarnings(as.integer(df_vcf[["POS"]]))

  # Check if there is NA in POS column after conversion
  invalid_pos_rows <- df_vcf[is.na(df_vcf[["POS"]]), ]

  # If line have no convertible POS value into int, print it
  if (nrow(invalid_pos_rows) > 0) {
    warning(paste0(
      "Warning: The following rows have non-integer values in the POS column: ",
      paste(invalid_pos_rows$row, collapse = ", ")
    ))
  }

  # Remove empty lines
  df_vcf <- df_vcf[rowSums(is.na(df_vcf) | df_vcf == "") != ncol(df_vcf), ]

  if (nrow(df_vcf) == 0) {
    stop("Error: The VCF file does not contain any mutation data.")
  }

  df_vcf
}



#' Parser TSV
#'
#' @param tsv_file Path to the TSV file.
#'
#' @return A cleaned data frame with required columns.
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
        na = character()
      )
    },
    error = function(e) {
      stop("Error: Failed to read the TSV file. Ensure it's properly formatted.")
    }
  )

  # Check if the read.table worked propely
  if (ncol(df_tsv) != 4) {
    stop(paste(
      "Error: The number of columns in the TSV file does not match the expected columns CHR POS REF ALT.",
      "Found columns:", paste(colnames(df_tsv), collapse = ", ")
    ))
  }

  # Convert POS column into int (Force to NA if str)
  df_tsv[["POS"]] <- suppressWarnings(as.integer(df_tsv[["POS"]]))

  # Check if there is NA in POS column after conversion
  invalid_pos_rows <- df_tsv[is.na(df_tsv[["POS"]]), ]

  # If line have no convertible POS value into int, print it
  if (nrow(invalid_pos_rows) > 0) {
    warning(paste0(
      "Warning: The following rows have non-integer values in the POS column: ",
      paste(invalid_pos_rows$row, collapse = ", ")
    ))
  }

  # Remove empty rows
  df_tsv <- df_tsv[rowSums(is.na(df_tsv) | df_tsv == "") != ncol(df_tsv), ]

  # Check if the TSV has data
  if (nrow(df_tsv) == 0) {
    stop("Error: The TSV file does not contain any mutation data.")
  }

  df_tsv
}


#' Parser chr:pos:ref:alt
#'
#' @param mut String chr:pos:ref:alt
#'
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
  if (!grepl("^[0-9]+$", parts[2])) {
    pos <- NA
  } else {
    pos <- suppressWarnings(as.integer(parts[2]))
  }

  ref <- ref <- parts[3]
  alt <- alt <- parts[4]

  # Return the info into a list
  mut_df <- data.frame(CHROM = chr, POS = pos, REF = ref, ALT = alt, stringsAsFactors = FALSE)

  # Check if the TSV has data
  if (nrow(mut_df) == 0) {
    stop("Error: The mutation string does not contain any mutation data.")
  }

  mut_df
}
