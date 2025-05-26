#' Get information about the repetition of an indel sequence
#'
#' Determines the sequence of an insertion or deletion (indel) from reference
#' and alternative alleles. It then scans the reference genome, starting
#' immediately after the indel's anchor base, to find how many times the
#' indel sequence unit repeats before a mismatch occurs or the chromosome ends.
#' @inheritParams process_fragmentomics
#' @param fasta_loaded An 'FaFile' object created by 'Rsamtools::FaFile()'
#'
#' @return A numeric value: the 1-based genomic position of the first nucleotide
#'   in the reference genome (downstream of 'pos') that does not match the
#'   pattern of the indel sequence.
#'   - If 'ref' and 'alt' have the same length (not an indel),
#'     it returns null and issues a warning.
#'   - If the derived indel sequence is empty, it returns null
#'     and issues a warning.
#'   - If the end of the chromosome is reached or a FASTA reading error occurs
#'     during comparison, it returns the genomic position where scanning stopped
#'
#' @importFrom Biostrings getSeq
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
#' @noRd
#' @noRd
get_repetition_seq_info <- function(chr, pos, ref, alt, fasta, mutation_type) {
  # Determine indel sequence based on user's VCF-like interpretation
  indel_seq <- ""

  if (mutation_type == "insertion") {
    indel_seq <- substring(alt, 2)
  } else if (mutation_type == "deletion") { # Deletion
    indel_seq <- substring(ref, 2)
  } else { # SNV or No change
    warning(paste0(
      "Input alleles (REF='", ref, "', ALT='", alt,
      "') do not represent an indel or mutation_type '", mutation_type, "' is not 'insertion' or 'deletion'. ",
      "This mutation cannot be processed as an indel by this function. ",
      "Returning NULL."
    ))
    return(NULL)
  }

  # Warning if the derived indel sequence is empty
  if (nchar(indel_seq) == 0) {
    warning(paste0(
      "Derived indel sequence is empty for REF='", ref, "', ALT='", alt,
      "', mutation_type='", mutation_type, "'. This might indicate an error in the indel representation ",
      "(e.g., alt/ref allele too short for substring(..., 2)). ",
      "This mutation cannot be processed as an indel. ",
      "Returning NULL."
    ))
    return(NULL)
  }

  # Initialize the current genomic position to the next base after the indel anchor
  current_genome_pos <- pos + 1
  indel_seq_len <- nchar(indel_seq)
  indel_char_idx <- 1 # 1-based index for substring

  # Max iterations to prevent potential infinite loops with unusual FASTA/indel inputs
  max_interations <- 100000
  iterations <- 0

  while (TRUE) {
    iterations <- iterations + 1
    if (iterations > max_interations) {
      warning(paste("Exceeded maximum iterations (", max_interations,
        ") at genomic position: ", current_genome_pos,
        ". Check for extremely long repeats or potential issues in FASTA/input.",
        sep = ""
      ))
      return(current_genome_pos)
    }

    # Retrieve the reference nucleotide from FASTA at current_genome_pos
    genome_base_seq_or_error <- tryCatch(
      {
        Biostrings::getSeq(
          x = fasta,
          param = GenomicRanges::GRanges(
            seqnames = chr,
            ranges = IRanges::IRanges(
              start = current_genome_pos,
              end = current_genome_pos
            )
          )
        )
      },
      error = function(e) {
        # Generate a more specific warning message
        warning_message <- paste(
          "Error trying to retrieve sequence from FASTA at",
          chr, ":", current_genome_pos, ".",
          "Details:", e$message
        )
        warning(warning_message)

        # Return a classed error object to differentiate error types
        # Check for common patterns in error messages for "chromosome not found"
        if (grepl("not found|non-existent sequence|no .faidx entry|failed to retrieve|record.*failed", e$message, ignore.case = TRUE)) {
          return(structure(list(message = e$message, pos = current_genome_pos), class = "ChromosomeNotFoundError"))
        } else {
          # For other errors like reading past end of an existing chromosome
          return(structure(list(message = e$message, pos = current_genome_pos), class = "SeqReadError"))
        }
      }
    )

    # Check the result of tryCatch
    if (inherits(genome_base_seq_or_error, "ChromosomeNotFoundError")) {
      return(NULL)
    } else if (inherits(genome_base_seq_or_error, "SeqReadError")) {
      return(current_genome_pos)
    }

    # Handle cases where getSeq() reads past the chromosome end without erroring
    # returning an empty DNAStringSet instead.
    if (length(genome_base_seq_or_error) == 0) {
      return(current_genome_pos)
    }

    # genome_base_seq_or_error is the actual DNAString sequence
    genome_base <- as.character(genome_base_seq_or_error)
    current_indel_base <- substring(indel_seq, indel_char_idx, indel_char_idx)

    if (genome_base != current_indel_base) {
      # Mismatch found
      return(current_genome_pos)
    }

    # Match: advance to the next base in genome and in indel sequence
    current_genome_pos <- current_genome_pos + 1
    indel_char_idx <- indel_char_idx + 1
    if (indel_char_idx > indel_seq_len) {
      indel_char_idx <- 1 # Cycle back to the start of the indel sequence
    }
  }
}
