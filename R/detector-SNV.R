# Project : ElsaBLab_fRagmentomics

#' Get base identity and sequencing quality from read and position of interest.
#'
#' @inheritParams build_fragments_info_table
#' @param r_pos numeric value representing the read mapping position
#' @param r_cigar character vector representing the read CIGAR
#' @param r_query character vector representing the read base sequence
#' @param r_qual character vector representing the read sequencing qualities
#' @return a named list with names `base` and `qual`.
#'
#' @importFrom stringr str_extract
#'
#' @noRd
get_base_qual_from_read <- function(pos, r_pos, r_cigar, r_query, r_qual){
  # table of operations-cursor movement
  Op <- c("M", "I", "D", "N", "S", "H", "P", "=", "X")
  Consumes_query <- c("yes", "yes", "no", "no", "yes", "no", "no", "yes", "yes")
  Consumes_reference <- c("yes", "no", "yes", "yes", "no", "no", "no", "yes", "yes")
  CIGAR <- data.frame(Op, Consumes_query, Consumes_reference)


  # cursors and index to be incremented
  cursor_reference <- r_pos
  index_query <- 1
  pos_is_readed <- FALSE

  while (nchar(r_cigar)>0 & cursor_reference<=pos){
    # get current cigar operation
    c_dg_str <- str_extract(r_cigar, "^[\\d]+")
    c_dg_num <- as.numeric(c_dg_str)
    c_op <- str_extract(r_cigar, "(?<=[\\d])[MIDNSHP=X]")
    r_cigar_left <- gsub(paste0("^", c_dg_str, c_op), "", r_cigar)


    # execute cigar operation to move cursors
    c_query <- CIGAR[CIGAR$Op==c_op, "Consumes_query"]
    c_reference <- CIGAR[CIGAR$Op==c_op, "Consumes_reference"]

    # If cursor reference is equal to the position of interest, then the cursor can advance to the position of interest
    if (cursor_reference == pos){
      pos_is_readed <- TRUE
    }

    # execute r_cigar operation along the reference
    if (c_reference=="yes"){
      cursor_reference <- cursor_reference + 1
    }

    #  execute r_cigar operation along the query
    if (c_query=="yes"){
      # we consume query and retriev new base and qual
      c_base <- substr(r_query, index_query, index_query)
      c_qual <- substr(r_qual, index_query, index_query)
      index_query <- index_query + 1
    } else {
      # we dont' consume query
      # reasons are: deletion (D), skipped base (N), hard clipping (H), padding (P)
      if (c_op=="D"){
        c_base <- "-"
        c_qual <- NA
      } else {
        c_base <- c_op
        c_qual <- NA
      }
    }

    # update cigar
    c_dg_num <- c_dg_num-1
    c_dg_str <- as.character(c_dg_num)
    if (c_dg_num==0){
      r_cigar <- r_cigar_left
    } else {
      r_cigar <- paste0(c_dg_str, c_op, r_cigar_left)
    }
  }

  # if the read does not cover the position of interest
  if (pos_is_readed == FALSE){
    c_base <- NA
    c_qual <- NA
  }

  c_indel <- NA

  list(base=c_base, qual=c_qual, indel=c_indel)
}