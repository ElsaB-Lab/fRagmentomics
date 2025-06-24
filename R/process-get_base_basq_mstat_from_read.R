#' Get read sequence, read base qualities, and read mutational status
#'
#' @inheritParams extract_fragment_features
#' @param read_stats A list of read-level statistics.
#'
#' @return A list containing "base", "basq", "mstat".
#'
#' @keywords internal
get_base_basq_mstat_from_read <- function(chr, pos, ref, alt, read_stats, fasta_fafile = NULL, fasta_seq = NULL,
                                          cigar_free_indel_match = FALSE) {
  # get index in the read sequence aligning with pos
  read_index_at_pos <- get_index_aligning_with_pos(pos, read_stats)

  if (read_index_at_pos == -1) {
    # in case the read does not cover the position of interest
    mstat <- NA
    base <- NA
    basq <- NA
  } else {
    # mutation status
    if (read_index_at_pos == -2 && nchar(ref) == nchar(alt)) {
      # in case the read contains a deletion at the position of interest and we are looking for a SNV or MNV,
      # then the read contains another event than the one we are looking for
      mstat <- "OTH (DEL)"
    } else if (read_index_at_pos == -2 && nchar(ref) != nchar(alt)) {
      mstat <- "OTH (DEL) for indel case - Need review"
    } else {
      # get mutation status of the read

      if (nchar(ref) == nchar(alt)) {
        # SNV or MNV
        mstat_small <- get_mutation_status_of_read(chr, pos, ref, alt, read_stats, read_index_at_pos, fasta_fafile,
          fasta_seq, cigar_free_indel_match,
          n_match_base_before = 0, n_match_base_after = 0
        )
        mstat_large <- get_mutation_status_of_read(chr, pos, ref, alt, read_stats, read_index_at_pos,
          fasta_fafile, fasta_seq, cigar_free_indel_match,
          n_match_base_before = 1, n_match_base_after = 1
        )
        if (mstat_small == "MUT") {
          if (mstat_large == "OTH") {
            mstat <- "MUT but potentially larger MNV"
          } else if (mstat_large == "MUT") {
            mstat <- "MUT"
          } else {
            stop(paste(
              "If the mutation status on the alt sequence is 'MUT', then the mutation status on extended",
              "sequence cannot be", paste0("'", mstat_large, "'"), "for",
              paste0(chr, ":", pos, ":", ref, ">", alt)
            ))
          }
        } else {
          mstat <- mstat_small
        }
      } else {
        # By VCF convention, the ref and alt sequences include one base before the actual mutated sequence
        # Additionally, we want to systematically include with one base after the actual mutated sequences.
        # If we cannot cover the base after, then we will call ambiguous if compatible with mutated ref.
        mstat <- get_mutation_status_of_read(chr, pos, ref, alt, read_stats, read_index_at_pos, fasta_fafile,
          fasta_seq, cigar_free_indel_match,
          n_match_base_before = 1, n_match_base_after = 1
        )
      }
    }

    # get bases to report. The idea is to report the bases aligning with pos, [inserted bases], pos+1, [inserted bases],
    # ..., pos+nchar(alt) and stop reporting bases once we reach nchar(alt)

    # initialize
    shift <- 0
    if (read_index_at_pos == -2) {
      base <- "-"
      basq <- ""
    } else {
      base <- substr(read_stats$SEQ, read_index_at_pos, read_index_at_pos)
      basq <- substr(read_stats$QUAL, read_index_at_pos, read_index_at_pos)
    }

    # iterate
    while (nchar(base) < nchar(alt)) {
      shift <- shift + 1
      info_at_pos <- get_base_basq_from_read_at_pos(pos + shift, pos, read_stats)
      max_new <- min(nchar(alt) - nchar(base), nchar(info_at_pos$base))
      base <- paste0(base, substr(info_at_pos$base, 1, max_new))
      basq <- paste0(basq, substr(info_at_pos$basq, 1, max_new))
      if (info_at_pos$base == "*") {
        break
      }
    }
  }

  list(base = base, basq = basq, mstat = mstat)
}


#' Get mutation status of read by comparing to wild-type and mutated reference
#'
#' @param read_seq Base sequence of the read
#' @param ref_seq_wt Base sequence of the wild-type ref
#' @param ref_seq_mut Base sequence of the mutated ref
#' @param compare_len_wt Number of bases to be compared in read vs wild-type ref comparison
#' @param compare_len_mut Number of bases to be compared in read vs mutated ref comparison
#' @return  a character describing the mutation status of the read
#'
#' @keywords internal
compare_read_to_ref_wt_and_mut <- function(read_seq, ref_seq_wt, ref_seq_mut, compare_len_wt, compare_len_mut) {
  # build the sequences to be compared
  ref_seq_wt_sub <- substr(ref_seq_wt, 1, compare_len_wt)
  ref_seq_mut_sub <- substr(ref_seq_mut, 1, compare_len_mut)
  read_seq_wt_sub <- substr(read_seq, 1, compare_len_wt)
  read_seq_mut_sub <- substr(read_seq, 1, compare_len_mut)

  # check if read compatible with wild-type ref
  match_ref_wt <- ref_seq_wt_sub == read_seq_wt_sub
  # check if read compatible with mutated ref
  match_ref_mut <- ref_seq_mut_sub == read_seq_mut_sub

  if (match_ref_wt && match_ref_mut) {
    return("AMB")
  } else if (match_ref_wt) {
    return("WT")
  } else if (match_ref_mut) {
    return("MUT")
  } else {
    return("OTH")
  }
}

#' Get the number of common first characters beetween two strings
#'
#'
#' @param str_a a string
#' @param str_b a string
#' @return An integer
#'
#' @keywords internal
get_number_of_common_first_char <- function(str_a, str_b) {
  i <- 0
  while (substr(str_a, i + 1, i + 1) == substr(str_b, i + 1, i + 1)) {
    i <- i + 1
  }

  i
}


#' Get the index of the element in the sequence aligning with the position of interest.
#'
#' Extracts base, quality, and indel informations based on the mutation type.
#'
#' @inheritParams get_base_basq_mstat_from_read
#' @return An integer representing the index of the nucleotide in sequence aligning with the position of interest. If
#' the read does not cover the position of interest, this integer is -1. If the read contains a deletion or skipping at
#' the position of interest, this integer is -2.
#'
#' @keywords internal
get_index_aligning_with_pos <- function(pos, read_stats) {
  # get POS, CIGAR, SEQ, and QUAL from the read
  read_pos <- read_stats$POS
  read_cigar <- read_stats$CIGAR

  # table of operations-cursor movement
  operations <- c("M", "I", "D", "N", "S", "H", "P", "=", "X")
  consumes_seq <- c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE)
  consumes_ref <- c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE)
  cigar_operations <- data.frame(operations, consumes_seq, consumes_ref)

  # Cursors to track our position along the reference and the read
  # ref_pos starts at the read's starting position (1-based)
  # read_idx starts at the first base of the query sequence (1-based)
  ref_pos <- read_pos
  read_idx <- 1

  while (nchar(read_cigar) > 0) {
    # get current cigar operation and number of associated bases
    read_cigar_nb_str <- str_extract(read_cigar, "^[\\d]+")
    read_cigar_nb <- as.numeric(read_cigar_nb_str)
    read_cigar_op <- str_extract(read_cigar, "(?<=[\\d])[MIDNSHP=X]")

    # get reduced cigar without the current operation
    read_cigar <- gsub(paste0("^", read_cigar_nb_str, read_cigar_op), "", read_cigar)

    # get properties of the current operation
    op_consumes <- cigar_operations[cigar_operations$operations == read_cigar_op, ]
    consumes_seq <- op_consumes$consumes_seq
    consumes_ref <- op_consumes$consumes_ref

    # Calculate how much this entire block moves the cursors
    ref_move <- if (consumes_ref) read_cigar_nb else 0
    seq_move <- if (consumes_seq) read_cigar_nb else 0

    # Check if the target position falls within the genomic block
    # covered by this operation. This only applies to ops that consume the reference.
    if (consumes_ref) {
      if (pos >= ref_pos && pos < (ref_pos + ref_move)) {
        # The position of interest is within this block.
        if (read_cigar_op %in% c("M", "=", "X")) {
          # The position is a match/mismatch. Find the exact index in the read.
          offset <- pos - ref_pos
          return(read_idx + offset)
        } else if (read_cigar_op %in% c("D", "N")) {
          # The position is a deletion or a skipped region in the read.
          # As per original logic, return -2.
          return(-2)
        }
      }
    }

    # If the position was not in the current block, update the cursors
    # to point to the beginning of the next block.
    ref_pos <- ref_pos + ref_move
    read_idx <- read_idx + seq_move
  }

  # if the loop finishes, the read did not cover the position of interest
  return(-1)
}


#' Get the element in the sequence and quality aligning with the position of interest.
#'
#' Extracts base, quality, and indel informations based on the mutation type.
#'
#' @inheritParams get_base_basq_mstat_from_read
#' @param pos_cur current position request
#' @param pos_pre previous position request
#' @return A list with the base and the quality of the read aligning with the position provided. If the read contains a
#' deletion, base is set to '-' and the quality to an empty string.
#'
#' @keywords internal
get_base_basq_from_read_at_pos <- function(pos_cur, pos_pre, read_stats) {
  read_index_at_pos_cur <- get_index_aligning_with_pos(pos_cur, read_stats)
  read_index_at_pos_pre <- get_index_aligning_with_pos(pos_pre, read_stats)

  if (read_index_at_pos_cur == -1) {
    return(list(base = "*", basq = "*"))
  } else if (read_index_at_pos_cur == -2) {
    return(list(base = "-", basq = ""))
  } else {
    base <- substr(read_stats$SEQ, read_index_at_pos_pre + 1, read_index_at_pos_cur)
    basq <- substr(read_stats$QUAL, read_index_at_pos_pre + 1, read_index_at_pos_cur)
    return(list(base = base, basq = basq))
  }
}
