# Project : ElsaBLab_fRagmentomics

#' @title bitvalues_from_bam_flag
#' @description Get bit values from bit names of BAM flag.
#'
#' @param flag an integer value
#' @param bitnames an integer vector
#' @return a matrix with column names matching bitnames and binary values representing bit values
#'
#' @keywords internal
bitvalues_from_bam_flag <- function(flag, bitnames) {
  FLAG_BITNAMES <- c(
      "isPaired",
      "isProperPair",
      "isUnmappedQuery",
      "hasUnmappedMate",
      "isMinusStrand",
      "isMateMinusStrand",
      "isFirstMateRead",
      "isSecondMateRead",
      "isSecondaryAlignment",
      "isNotPassingQualityControls",
      "isDuplicate",
      "isSupplementaryAlignment"
  )


  bitpos <- match(bitnames, FLAG_BITNAMES)
  invalid_bitnames_idx <- which(is.na(bitpos))
  if (length(invalid_bitnames_idx) != 0L) {
      in1string <- paste0(bitnames[invalid_bitnames_idx], collapse=", ")
      stop("invalid bitname(s): ", in1string)
  }

  ans <- S4Vectors:::explodeIntBits(flag, bitpos=bitpos)
  dimnames(ans) <- list(names(flag), bitnames)
  ans
}

#' @title parse_cigar
#' @description Function to parse the CIGAR string
#'
#' @param cigar a character vector
#' @return a df with length and type of the CIGAR string
#'
#' @keywords internal
parse_cigar <- function(cigar) {
  # Extract the operations and lengths from the CIGAR string
  matches <- gregexpr("([0-9]+)([MIDNSHP=X])", cigar, perl=TRUE)
  ops_str <- regmatches(cigar, matches)[[1]]

  lengths <- as.numeric(gsub("([0-9]+)([MIDNSHP=X])", "\\1", ops_str))
  types <- gsub("([0-9]+)([MIDNSHP=X])", "\\2", ops_str)

  data.frame(length = lengths, type = types, stringsAsFactors = FALSE)
}


#' @title Check sequence repetition
#' @description Check if the indel is in a repeted sequence
#' 
#' @param indel_rep an integer value
#' @param is_length_ok an booleen value
#' @param diff_mutation_bam an integer value
#' @param deletion_length an integer value
#' @return 0 or 1 if the indel is the one we're looking for
#'
#' @keywords internal
check_seq_rep <- function(indel_rep, is_length_ok, diff_mutation_bam, deletion_length) {
  if (!is_length_ok) {
    return(0)
  }
  if (diff_mutation_bam == 0) {
    return(1)
  }
  if ((diff_mutation_bam %% deletion_length) == 0) {
    factor <- abs(diff_mutation_bam) / deletion_length
    return(as.integer(factor <= indel_rep))
  }
  return(0)
}


#' @title Define the number of times the insertion is repeted
#' @description Return an integer value which represents the number of the sequence insertion that is repeted in the sequence 
#'
#' @param alt_len an integer value representing the length of the insertion 
#' @param current_pos an integer value
#' @param r_pos numeric value representing the read mapping position
#' @param cigar character vector representing the read CIGAR
#' @param r_query character vector representing the read base sequence
#' @return an integer value which represents the number of the sequence insertion that is repeted in the sequence 
#'
#' @keywords internal
define_ins_rep <- function(alt_len, current_pos, r_pos, cigar, r_query) {
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
  insertion_seq <- substr(r_query, pos_insertion + 1, pos_insertion + alt_len)

  # Move leftward in steps of alt_len to find the leftmost starting position of the insertion sequence
  pos_left <- pos_insertion + 1
  while (pos_left - alt_len >= 1) {
    # Extract the sequence of length alt_len to the left
    candidate <- substr(r_query, pos_left - alt_len, pos_left - 1)
      
    if (candidate == insertion_seq) {
      # If the extracted sequence matches the insertion sequence, move left
      pos_left <- pos_left - alt_len
    } else {
      # Stop if the sequence doesn't match
      break
    }
  }

  # Now, pos_left is the leftmost position where the repetition starts.
  # Count how many times the insertion sequence repeats starting from pos_left
  ins_rep <- 0
  pos_check <- pos_left

  while (pos_check + alt_len - 1 <= nchar(r_query)) {
    # Extract the sequence of length alt_len to the right
    candidate <- substr(r_query, pos_check, pos_check + alt_len - 1)
      
    if (candidate == insertion_seq) {
      # If the extracted sequence matches, increment the count and move right
      ins_rep <- ins_rep + 1
      pos_check <- pos_check + alt_len
    } else {
      # Stop if the sequence doesn't match
      break
    }
  }
  return(ins_rep)
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

  # Keep the sequence from the beggining of the insertion 
  insertion_qual <- substr(r_qual, pos_insertion + 1, pos_insertion + alt_len)
  
  return(insertion_qual)
}