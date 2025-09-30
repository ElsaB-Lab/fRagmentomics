#' Extract allele, quality, and mutation status from a read
#'
#' @description This function serves as a primary parser for a single read at a variant locus. It determines if the
#' read covers the variant, extracts the observed sequence and its base qualities, and assigns a detailed mutation status.
#'
#' @inheritParams extract_fragment_features
#' @param read_stats A list of read-level statistics.
#'
#' @return A list containing "base", "basq", "mstat".
#'
#' @keywords internal
get_base_basq_mstat_from_read <- function(chr, pos, ref, alt, read_stats, fasta_fafile = NULL, fasta_seq = NULL) {
  # get index in the read sequence aligning with pos
  read_index_at_pos <- get_index_aligning_with_pos(pos, read_stats)

  if (read_index_at_pos %in% c(-1, -0.5)) {
    # in case the read does not cover the position of interest
    mstat <- NA
    base <- NA
    basq <- NA
  } else {
    # mutation status
    if (read_index_at_pos == -2) {
      # in case the read contains a deletion at the position of interest then the read contains another event than the
      # one we are looking for
      mstat <- "OTH"
    } else {
      # get mutation status of the read

      if (nchar(ref) == nchar(alt)) {
        # SNV or MNV
        mstat_small <- get_mutation_status_of_read(chr, pos, ref, alt, read_stats, read_index_at_pos, fasta_fafile,
          fasta_seq,
          n_match_base_before = 0, n_match_base_after = 0
        )
        mstat_large <- get_mutation_status_of_read(chr, pos, ref, alt, read_stats, read_index_at_pos,
          fasta_fafile, fasta_seq,
          n_match_base_before = 1, n_match_base_after = 1
        )
        if (mstat_small == "MUT") {
          if (mstat_large == "OTH") {
            mstat <- "MUT but potentially OTH"
          } else if (mstat_large == "MUT") {
            mstat <- "MUT"
          } else {
            stop(sprintf(
              "If the mutation status on the alt sequence is 'MUT', then the mutation status on extended sequence cannot be '%s' for %s:%d:%s>%s",
              mstat_large,
              chr,
              pos,
              ref,
              alt
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
          fasta_seq,
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
      info_at_pos <- get_base_basq_from_read_at_pos(pos + shift, pos + shift - 1, read_stats)
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


#' Extract sequence and quality between two read indices
#'
#' @description Extracts the read sequence and base qualities corresponding to the region between two read indices,
#' which are derived from genomic positions.
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
  } else if (read_index_at_pos_cur == -0.5) {
    base <- substr(read_stats$SEQ, read_index_at_pos_pre + 1, nchar(read_stats$SEQ))
    basq <- substr(read_stats$QUAL, read_index_at_pos_pre + 1, nchar(read_stats$QUAL))
    return(list(base = base, basq = basq))
  } else {
    base <- substr(read_stats$SEQ, read_index_at_pos_pre + 1, read_index_at_pos_cur)
    basq <- substr(read_stats$QUAL, read_index_at_pos_pre + 1, read_index_at_pos_cur)
    return(list(base = base, basq = basq))
  }
}
