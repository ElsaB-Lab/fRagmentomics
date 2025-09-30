#' Parse a CIGAR string into its component operations
#'
#' @description
#' This internal helper function deconstructs a CIGAR string into a structured
#' format, separating the operation type (e.g., M, I, D) from its length.
#'
#' @param cigar a character vector
#'
#' @return a df with length and type of the CIGAR string. One row per operation.
#'
#' @noRd
parse_cigar <- function(cigar) {
  # Extract the operations and lengths from the CIGAR string
  matches <- gregexpr("([0-9]+)([MIDNSHP=X])", cigar, perl = TRUE)
  ops_str <- regmatches(cigar, matches)[[1]]

  lengths <- as.numeric(gsub("([0-9]+)([MIDNSHP=X])", "\\1", ops_str))
  types <- gsub("([0-9]+)([MIDNSHP=X])", "\\2", ops_str)

  data.frame(length = lengths, type = types, stringsAsFactors = FALSE)
}

#' Calculate read length excluding trailing soft clips
#'
#' @description
#' This function computes the effective length of a read by subtracting the
#' number of bases that are soft-clipped at its 3' end, which is a common
#' step in variant analysis to avoid alignment artifacts.
#'
#' @param cigar A CIGAR string
#' @param seq The corresponding DNA sequence
#'
#' @return The length of the sequence after subtracting the number of bases
#'         that are soft-clipped at the end of the read.
#'
#' @noRd
calculate_len_without_end_softclip <- function(cigar, seq) {
  read_len <- nchar(seq)
  softclip_len <- 0

  # Check if the CIGAR string ends with a soft clip and capture the number of bases
  if (grepl("\\d+S$", cigar)) {
    # Extract the numeric value of the trailing soft clip
    softclip_len <- as.numeric(sub(".*?(\\d+)S$", "\\1", cigar))
  }

  # Return the total length minus the trailing soft clip length
  return(read_len - softclip_len)
}

#' Extract insertion and deletion positions from a read
#'
#' @description
#' Parses a read's CIGAR string and starting position to determine the precise
#' genomic coordinates of any insertions or deletions.
#'
#' @param read_stats A list containing read infos.
#' @return A list containing two vectors: deletions and insertions.
#'
#' @noRd
get_pos_indels_from_read <- function(read_stats) {
  list_pos_del <- c()
  list_pos_ins <- c()

  # Parse the CIGAR string
  cigar_ops <- parse_cigar(read_stats$CIGAR)

  # Initialize the current genomic position
  current_ref_pos <- read_stats$POS

  # Iterate over each CIGAR operation
  for (i in seq_len(nrow(cigar_ops))) {
    op_len <- cigar_ops$length[i]
    op_type <- cigar_ops$type[i]

    # Handle operations that consume the reference
    if (op_type %in% c("M", "=", "X")) {
      # Match or Mismatch: advances the position on the reference
      current_ref_pos <- current_ref_pos + op_len
    } else if (op_type == "D") {
      # Deletion: does not advance on the read, but does on the reference
      # The positions of the deleted bases
      deleted_positions <- seq(from = current_ref_pos, length.out = op_len)
      list_pos_del <- c(list_pos_del, deleted_positions)

      # Update the position on the reference
      current_ref_pos <- current_ref_pos + op_len
    } else if (op_type == "I") {
      # Insertion: advances on the read, but not on the reference
      # The genomic position is the one just before the insertion
      pos_before_insertion <- current_ref_pos - 1

      # Create the custom positions for each inserted base
      inserted_positions <- as.numeric(paste0(pos_before_insertion, ".", seq_len(op_len)))
      list_pos_ins <- c(list_pos_ins, inserted_positions)
      # current_ref_pos does not change
    } else if (op_type == "N") {
      # Skipped region: advances the position on the reference
      current_ref_pos <- current_ref_pos + op_len
    }
    # 'S' (soft clip), 'H' (hard clip), and 'P' (padding) operations
    # do not consume the reference, so we do nothing.
  }

  return(list(deletions = list_pos_del, insertions = list_pos_ins))
}

#' Create an empty data row for a QC-failed fragment
#'
#' @description
#' Generates a standardized single-row 'data.frame' for a DNA fragment that
#' failed quality control. All data fields are populated with 'NA', ensuring
#' consistent output structure even for excluded fragments.
#'
#' @inheritParams extract_fragment_features
#' @param fragment_qc a character vector corresponding to QC failure of the fragment
#'
#' @return a df with length and type of the CIGAR string. One row per operation.
#'
#' @noRd
create_empty_fragment_row <- function(
    chr, pos, ref, alt, input_mutation_info, fragment_name, fragment_qc,
    sample_id, report_tlen, report_5p_3p_bases_fragment, report_softclip) {
  final_row_fragment <- list(
    Chromosome             = chr,
    Position               = pos,
    Ref                    = ref,
    Alt                    = alt,
    Input_Mutation         = input_mutation_info,
    Fragment_Id            = fragment_name,
    Fragment_QC            = fragment_qc,
    Fragment_Status_Simple = NA_character_,
    Fragment_Status_Detail = NA_character_,
    Fragment_Size          = NA_integer_,
    Read_5p_Status         = NA_character_,
    Read_3p_Status         = NA_character_,
    FLAG_5p                = NA_integer_,
    FLAG_3p                = NA_integer_,
    MAPQ_5p                = NA_integer_,
    MAPQ_3p                = NA_integer_,
    BASE_5p                = NA_character_,
    BASE_3p                = NA_character_,
    BASQ_5p                = NA_character_,
    BASQ_3p                = NA_character_,
    CIGAR_5p               = NA_character_,
    CIGAR_3p               = NA_character_,
    POS_5p                 = NA_integer_,
    POS_3p                 = NA_integer_
  )

  if (!is.na(sample_id)) {
    final_row_fragment$Sample_Id <- sample_id
  }

  if (report_tlen) {
    final_row_fragment$TLEN <- NA_integer_
  }

  if (report_5p_3p_bases_fragment != 0) {
    final_row_fragment$Fragment_Bases_5p <- NA_character_
    final_row_fragment$Fragment_Bases_3p <- NA_character_
    final_row_fragment$Fragment_Basqs_5p <- NA_character_
    final_row_fragment$Fragment_Basqs_3p <- NA_character_
  }

  if (report_softclip) {
    final_row_fragment$Nb_Fragment_Bases_Softclip_5p <- NA_integer_
    final_row_fragment$Nb_Fragment_Bases_Softclip_3p <- NA_integer_
  }

  # Add VAF column
  final_row_fragment$VAF <- NA_real_

  result_df <- as.data.frame(final_row_fragment)
  return(result_df)
}

#' Save a plot to disk if 'output_path' is provided
#'
#' @description
#' Internal helper to save a ggplot object using 'ggplot2::ggsave()'.
#' Validates 'output_path', requires an existing directory, and merges
#' user params with sensible defaults.
#'
#' @param final_plot a ggplot object to save.
#' @param output_path single string; full file path. If 'NULL'/empty, no save occurs.
#' @param ggsave_params named list of extra args for 'ggplot2::ggsave()' (overrides defaults).
#' @param defaults named list of default 'ggsave()' args (width, height, units, dpi, bg).
#'
#' @return returns 'final_plot' (invisibly) if not saved; 'NULL' (invisibly) after saving.
#'
#' @keywords internal
#' @noRd
#' @keywords internal
save_plot_if_needed <- function(final_plot, output_path, ggsave_params = list(),
                                 defaults = list(width = 8, height = 6, units = "in", dpi = 300, bg = "white")) {
  # Check if the file already exists and warn the user if it will be overwritten
  if (file.exists(output_path)) {
    message(sprintf("File '%s' already exists and will be overwritten.", output_path))
  } else {
    # Check if the file parent folder already exists
    output_parent <- dirname(output_path)
    if (!dir.exists(output_parent)){
      message(sprintf("Folder '%s' does not exist and will be created", output_parent))
      dir.create(output_parent, showWarnings=FALSE, recursive=TRUE)
    }
  }

  # Merge params (user overrides)
  final_params <- utils::modifyList(defaults, ggsave_params)
  args <- c(list(plot = final_plot, filename = output_path), final_params)

  message(sprintf("Saving plot to: %s", output_path))

  tryCatch(
    do.call(ggplot2::ggsave, args),
    error = function(e) stop(sprintf("Failed to save plot to '%s': %s", output_path, e$message), call. = FALSE)
  )

  invisible(NULL)
}

