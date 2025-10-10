#' Get mutation status of read
#'
#' @description
#' The algorithm to determine if a read supports the mutation of interest is quite intuitive but requires careful
#' considerations of insertions and deletions in repeated sequences to resolve ambiguities whenever possible. The core
#' idea is to apply the mutation on the reference sequence and compare the mutated reference sequence to the read
#' sequence starting from the position of the mutation of interest. The sequence of the read should be compared to the
#' sequence of the reference sequence and of the mutated reference sequence for just enough bases to assign mutational
#' status without ambiguity.
#'
#' @inheritParams get_base_basq_mstat_from_read
#' @param read_index_at_pos An integer representing the index of the nucleotide in sequence aligning with the position of interest.
#' @param n_match_base_before Number of bases to be matched before the alt allele in the sequences comparison
#' @param n_match_base_after Number of bases to be matched after the last alt allele in the sequences comparison
#'
#' @return
#' A character string indicating the mutational status of the read. Possible values include:
#' \itemize{
#'  \item "WT": Wild-Type. The read matches the reference sequence.
#'  \item "MUT": Mutant. The read matches the alternate allele.
#'  \item "AMB": Ambiguous. The read is too short or the context is too complex to definitively assign a status.
#'  \item "OTH": Other. An alteration is found, but it is not the one of interest.
#'  \item For INDELs, the status may be combined with a descriptive message (e.g., "AMB by cigar-free search...").
#' }
#'
#' @keywords internal
get_mutation_status_of_read <- function(chr, pos, ref, alt, read_stats, read_index_at_pos, fasta_fafile = NULL,
                                        fasta_seq = NULL, n_match_base_before = 1,
                                        n_match_base_after = 1) {
  ref_len <- nchar(ref)
  alt_len <- nchar(alt)

  read_seq_len <- nchar(read_stats$SEQ)

  # first identify the minimum number of bases that should be considered in the comparison
  # this number should be at the very least the size of the one base before + size of the mutated sequence  + one base
  # after
  # in case of an insertion or a deletion, additional bases should be included in the comparison if the reference
  # has a repeat of the mutated sequence at the locus considered
  if (ref_len == alt_len) {
    # SNV or MNV

    # If the read sequence does not allow to compare with base(s) before the alteration, or with base(s) after the
    # alteration, then we won't add these bases into the comparison
    if (read_index_at_pos == 1) {
      n_match_base_before <- 0
    }
    if ((read_index_at_pos + alt_len - 1) >= read_seq_len) {
      n_match_base_after <- 0
    }

    # we fetch one base before, the size of the allele, and one base after
    fetch_len_ref <- n_match_base_before + alt_len + n_match_base_after
    fetch_start_ref <- pos - n_match_base_before
    fetch_end_ref <- fetch_start_ref + fetch_len_ref - 1

    # get sequences of wild-type and mutated ref
    ref_seq_wt <- get_seq_from_fasta(chr, fetch_start_ref, fetch_end_ref, fasta_fafile, fasta_seq)
    ref_seq_mut <- paste0(
      substr(ref_seq_wt, 1, n_match_base_before), alt,
      substr(ref_seq_wt, nchar(ref_seq_wt) - n_match_base_after + 1, nchar(ref_seq_wt))
    )

    # if we need more bases than available, cut the sequences to only the maximum comparison possible
    # for instance
    # REF: ATTCGAGTAT
    # At position 9, AT > CC
    # READ_SEQ: AGTCC
    # Then,
    #   read_index_at_pos: 4
    #   fetch_start_read: 3
    #   read_seq_len: 5
    #   fetch_len: 3 (as we can match one base before but can't match any base after)
    # We want to fetch the bases 3,4,5
    fetch_start_read <- read_index_at_pos - n_match_base_before
    fetch_len_read <- min(fetch_len_ref, read_seq_len - fetch_start_read + 1)
    fetch_end_read <- fetch_start_read + fetch_len_read - 1
    read_seq <- substr(read_stats$SEQ, fetch_start_read, fetch_end_read)

    # run comparison on maximum size possible
    compare_len_wt <- n_match_base_before + alt_len + n_match_base_after
    compare_len_mut <- n_match_base_before + alt_len + n_match_base_after

    # if we cannot cover completely the allele with the read, we may have an ambiguity
    incomplete_comparison_mut <- (fetch_len_read < compare_len_mut)
    compare_len_wt <- min(compare_len_wt, fetch_len_read)
    compare_len_mut <- min(compare_len_mut, fetch_len_read)

    # determine the read mutation status
    status <- compare_read_to_ref_wt_and_mut(read_seq, ref_seq_wt, ref_seq_mut, compare_len_wt, compare_len_mut)

    if (status == "MUT" && incomplete_comparison_mut) {
      # An incomplete comparison compatible with the mutated allele and not compatible with the reference allele
      # will be labelled "AMB" because we don't observe completely the alternative allele
      return("AMB")
    } else {
      return(status)
    }
  } else {
    # INS or DEL

    # identify mutated sequence
    inserted_seq <- if (alt_len > ref_len) substr(alt, 2, alt_len) else ""
    deleted_seq <- if (ref_len > alt_len) substr(ref, 2, ref_len) else ""
    motif <- if (inserted_seq != "") inserted_seq else if (deleted_seq != "") deleted_seq
    motif_len <- nchar(motif)

    # Fetch the maximum number of bases that may be included in the comparison
    fetch_len_ref <- n_match_base_before - 1 + read_seq_len - (read_index_at_pos - 1) + motif_len + n_match_base_after
    fetch_start_ref <- pos - (n_match_base_before - 1)
    fetch_end_ref <- fetch_start_ref + fetch_len_ref - 1
    ref_seq_wt <- get_seq_from_fasta(chr, fetch_start_ref, fetch_end_ref, fasta_fafile, fasta_seq)

    # INS example
    #   REF: TGAGAT
    #   MUT: T > TGA
    #   REF_AFTER_ONE_MOTIF: GAT
    #
    # DEL example
    #   REF: TGAGAGAT
    #   MUT: TGA > T
    #   REF_AFTER_ONE_MOTIF: GAGAT

    # if the sequence contains the searched motif at least once
    regex_motif <- paste0("^", if (alt_len > ref_len) alt else ref)

    if (grepl(regex_motif, substr(ref_seq_wt, n_match_base_before, nchar(ref_seq_wt)))) {
      # Because we are supposed to find a deletion motif and
      repeat_count <- 1
      ref_seq_after_one_motif <- substr(ref_seq_wt, n_match_base_before + motif_len + 1, nchar(ref_seq_wt))

      # identify the number of repeats within the region of the maximum comparison size
      while (substr(ref_seq_after_one_motif, (repeat_count - 1) * motif_len + 1, repeat_count * motif_len) == motif) {
        repeat_count <- repeat_count + 1
      }

      # Identify how many additional bases need to be included in the comparison
      # This number is the number of bases of the motif that match with the reference sequence  after the last repeat
      # of the motif
      # As examples,
      #   - If REF is A A T C A T C A G T and we are looking for the insertion A > AATC
      #   Here the motif "ATC" is repeated twice. The sequence after the last repeat is AGT. This sequence shares one
      #   common base with the motif, i.e "A"
      #
      #   - If REF is A A T C A T C A T T and we are looking for the insertion A > AATC
      #   Here the motif "ATC" is repeated twice. The sequence after the last repeat is ATT. This sequence shares two
      #   common bases with the motif, i.e "AT"
      ref_seq_after_last_motif <- substr(
        ref_seq_after_one_motif, (repeat_count - 1) * motif_len + n_match_base_after,
        repeat_count * motif_len + n_match_base_after - 1
      )
      n_bases_shared_with_motif <- get_number_of_common_first_char(ref_seq_after_last_motif, motif)
    } else {
      repeat_count <- 0
      n_bases_shared_with_motif <- 0
    }

    # Compute the size of the comparisons to wild-type ref and mutated ref
    if (alt_len > ref_len) {
      type <- "I"
      compare_len_wt <- n_match_base_before + repeat_count * motif_len + n_bases_shared_with_motif + n_match_base_after
      compare_len_mut <- n_match_base_before + (repeat_count + 1) * motif_len + n_bases_shared_with_motif + n_match_base_after
    } else {
      type <- "D"
      compare_len_wt <- n_match_base_before + (repeat_count - 1) * motif_len + n_bases_shared_with_motif + n_match_base_after
      compare_len_mut <- n_match_base_before + (repeat_count - 1) * motif_len + n_bases_shared_with_motif + n_match_base_after
    }

    # Build mutated sequence by inserting or deleting the correct bases
    # INS example
    #   REF: TGAGAT
    #   MUT: T > TGA
    #   REF_SEQ_MUT: T+GA+GAGAT
    #
    # DEL example
    #   REF: TGAGAGAT
    #   MUT: TGA > T
    #   REF_SEQ_MUT: T+""+GAGAT
    ref_seq_mut <- paste0(
      substr(ref_seq_wt, 1, n_match_base_before),
      substr(alt, n_match_base_before + 1, alt_len),
      substr(ref_seq_wt, ref_len + 1, nchar(ref_seq_wt))
    )

    # Fetch maximum read sequence available starting from the position needed to cover n_match_base_before
    fetch_start_read <- read_index_at_pos - (n_match_base_before - 1)
    fetch_len_read <- read_seq_len - fetch_start_read + 1
    fetch_end_read <- fetch_start_read + fetch_len_read - 1
    read_seq <- substr(read_stats$SEQ, fetch_start_read, fetch_end_read)

    # if we cannot cover completely the allele with the read, we may have an ambiguity
    incomplete_comparison_mut <- (fetch_len_read < compare_len_mut)
    compare_len_wt <- min(compare_len_wt, fetch_len_read)
    compare_len_mut <- min(compare_len_mut, fetch_len_read)

    # determine the read mutation status
    status <- compare_read_to_ref_wt_and_mut(read_seq, ref_seq_wt, ref_seq_mut, compare_len_wt, compare_len_mut)
    indel_search <- search_for_indel_in_cigar(pos, ref, alt, read_stats, type)
    indel_found_in_cigar <- indel_search[[1]]
    other_found_in_cigar <- indel_search[[2]]

    if (!incomplete_comparison_mut) {
      # -------------------- complete_comparison --------------------
      if (indel_found_in_cigar) {
        switch(status,
          "MUT" = "MUT",
          "WT"  = "MUT by CIGAR but potentially WT",
          "AMB" = "IMPOSSIBLE",
          "OTH" = "MUT by CIGAR but potentially OTH"
        )
      } else if (other_found_in_cigar) {
        switch(status,
          "MUT" = "OTH by CIGAR but potentially MUT",
          "WT"  = "OTH by CIGAR but potentially WT",
          "AMB" = "IMPOSSIBLE",
          "OTH" = "OTH"
        )
      } else { # mutation not found by CIGAR
        switch(status,
          "MUT" = "WT by CIGAR but potentially MUT",
          "WT"  = "WT",
          "AMB" = "IMPOSSIBLE",
          "OTH" = "OTH"
        )
      }
    } else {
      # -------------------- incomplete_comparison --------------------
      if (indel_found_in_cigar) {
        switch(status,
          "MUT" = "MUT by CIGAR but AMB",
          "WT"  = "MUT by CIGAR but potentially WT",
          "AMB" = "MUT by CIGAR but AMB",
          "OTH" = "MUT by CIGAR but potentially OTH"
        )
      } else if (other_found_in_cigar) {
        switch(status,
          "MUT" = "OTH by CIGAR but potentially MUT",
          "WT"  = "OTH by CIGAR but potentially WT",
          "AMB" = "OTH by CIGAR but AMB",
          "OTH" = "OTH"
        )
      } else { # mutation not found by CIGAR
        switch(status,
          "MUT" = "WT by CIGAR but potentially MUT",
          "WT"  = "WT",
          "AMB" = "AMB",
          "OTH" = "OTH"
        )
      }
    }
  }
}


#' Count common leading characters between two strings
#'
#' @description
#' Calculates the length of the common prefix for two given character strings.
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
