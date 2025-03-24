# Project : ElsaBLab_fRagmentomics

#' sanity_check_ins
#'
#' @param alt_len an integer value representing the length of the insertion 
#' @param current_pos an integer value
#' @param r_pos numeric value representing the read mapping position
#' @param cigar character vector representing the read CIGAR
#' @param r_query character vector representing the read base sequence
#' @return a character vector representing the insertion sequence 
#'
#' @keywords internal
get_seq_ins <- function(alt_len, current_pos, r_pos, cigar, r_query) {
  # Extraction of the sum of the cigar with the S
  matches <- gregexpr("^[0-9]+S", cigar)  # The ^ means start of the chain 
  extracted <- regmatches(cigar, matches)[[1]]

  # Put 0 if no soft clipping
  if (length(extracted) == 0) {
    start_soft_clip <- 0  
  } else {
    start_soft_clip <- as.numeric(gsub("S", "", extracted))
  }
  
  pos_insertion <- start_soft_clip + current_pos - r_pos

  # Keep the sequence from the beggining of the insertion 
  insertion_seq <- substr(r_query, pos_insertion, pos_insertion + alt_len)

  return(insertion_seq)
}

#' @title Find quality of the base
#' @description Find quality of the base
#'
#' @param alt_len an integer value representing the length of the insertion 
#' @param current_pos an integer value representing the current position of the insertion before the cigar operations
#' @param r_pos numeric value representing the read mapping position
#' @param cigar character vector representing the read CIGAR
#' @param r_qual character vector representing the read sequencing qualities
#' @return 0 or 1 if the indel is the one we're looking for
#'
#' @keywords internal
find_c_qual <- function(alt_len, current_pos, r_pos, r_cigar, r_qual) {
  # Extraction of the sum of the cigar with the S
  matches <- gregexpr("^[0-9]+S", r_cigar)  # The ^ means start of the chain 
  extracted <- regmatches(r_cigar, matches)[[1]]

  # Put 0 if no soft clipping
  if (length(extracted) == 0) {
    start_soft_clip <- 0  
  } else {
    start_soft_clip <- as.numeric(gsub("S", "", extracted))
  }

  pos_insertion <- start_soft_clip + current_pos - r_pos 

  # Keep the sequence from the nucleotide before the insertion (like in ref)
  insertion_qual <- substr(r_qual, pos_insertion, pos_insertion + alt_len)
  
  return(insertion_qual)
}

#' Get the informations about the presence of a insertion in the read.
#'
#' @inheritParams fRagmentomics
#' @param pos numeric value representing the position of interest
#' @param alt character vector representing the insertion sequence
#' @param r_pos numeric value representing the read mapping position
#' @param r_cigar character vector representing the read CIGAR
#' @param r_query character vector representing the read base sequence
#' @param r_qual character vector representing the read sequencing qualities
#' @return a named list with names `base`,`qual`.
#'
#' @noRd
get_insertion <- function(pos, alt, r_pos, r_cigar, r_query, r_qual) {
  c_base <- "no_insertion_detected"
  c_qual <- "no_insertion_detected"
  pos_is_readed <- FALSE

  # Parse the CIGAR string
  ops <- parse_cigar(r_cigar)

  # Iterating over CIGAR operations
  current_pos <- r_pos  # Current reference position, aligned with the read
  for (i in seq_len(nrow(ops))) {
    op_len <- ops$length[i]
    op_type <- ops$type[i]

    if (op_type %in% c("M", "D", "N", "=", "X")) {
      # M: match or mismatch, D: deletion in the read, N: skip in the reference, =: perfect match, X: mismatch
      current_pos <- current_pos + op_len

    } else if (op_type == "I") {
      # I: insertion in the read (addition in the read compared to the reference)
      # Compare this insertion with the reference
      diff_mutation_bam <- pos - (current_pos - 1)

      # Sanity check of the insertion sequence
      ins_seq <- get_seq_ins(op_len, current_pos, r_pos, r_cigar, r_query)

      if ((ins_seq == alt) && (diff_mutation_bam == 0)){
        c_base <- "insertion_detected"
        c_qual <- find_c_qual(op_len, current_pos, r_pos, r_cigar, r_qual)
        pos_is_readed <- TRUE
        break
      }

    } else {
      # S, H, P: Soft clip, Hard clip, Pad - do not advance the reference position
      # No action is taken
    }
  }

  # Check if read covers the position of interest
  if ((current_pos - 1) > pos) {
    pos_is_readed <- TRUE
  }

  # if the read does not cover the position of interest
  if (pos_is_readed == FALSE){
    c_base <- NA
    c_qual <- NA
  }

  return(list(base = c_base, qual = c_qual))
}