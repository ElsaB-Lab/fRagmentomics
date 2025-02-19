# Project : cfDNA mutation origin


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


#' Get the informations about the presence of a deletion in the read.
#'
#' @inheritParams build_fragments_info_table*
#' @param chr numeric value representing the chromosome
#' @param pos numeric value representing the position of interest
#' @param ref character vector representing the deletion sequence
#' @param del_info character vector sous la forme "$pos_final,$rep_del"
#' @param r_pos numeric value representing the read mapping position
#' @param r_cigar character vector representing the read CIGAR
#' @return a named list with names `base`,`qual` and `indel`.
#'
#' @importFrom stringr str_extract utils
#'
#' @noRd
get_deletion <- function(chr, pos, ref, del_info, r_pos, r_cigar) {
  # Split the del_info parameter into final_pos and del_rep
  vars <- strsplit(del_info, ",")[[1]]
  final_pos <- suppressWarnings(as.numeric(vars[1]))
  del_rep <- as.numeric(vars[2])

  # Initialize variables
  c_base <- NA
  c_qual <- NA
  c_indel <- "No_support_del"   # by default, no deletion
  ref_len <- nchar(ref)

  # Check if final_pos is not numeric
  if (is.na(final_pos)) {
    c_indel <- "no_del_found_in_ref_genome"
    return(list(base = c_base, qual = c_qual, indel = c_indel))

  } else {
    # Parse the CIGAR string
    ops <- parse_cigar(r_cigar)

    # Iterating over CIGAR operations
    current_pos <- r_pos  # Current reference position, aligned with the read
    for (i in seq_len(nrow(ops))) {
      op_len <- ops$length[i]
      op_type <- ops$type[i]

      if (op_type == "M") {
        # M: match or mismatch (but advances the reference position)
        current_pos <- current_pos + op_len
      } else if (op_type == "I") {
        # I: insertion in the read (does not advance the reference position, only advances in the read)
        # Therefore, we do not modify current_pos
        # (current_pos remains the reference position)
      } else if (op_type == "D") {
        # D: deletion in the read (missing in the read compared to the reference)
        # Compare this deletion to the one in the reference
        diff_mutation_bam <- final_pos - (current_pos - 1)
        is_length_ok <- (ref_len == op_len)

        c_indel <- check_seq_rep(del_rep, is_length_ok, diff_mutation_bam, op_len)

        # Advance the reference position after the deletion
        current_pos <- current_pos + op_len

        if (c_indel == 1) {
          # A corresponding deletion was found
          break
        }
      } else if (op_type %in% c("N", "=", "X")) {
        # N: skip in the reference, =: perfect match, X: mismatch
        # All advance the reference position
        current_pos <- current_pos + op_len
      } else {
        # S, H, P: Soft clip, Hard clip, Pad - do not advance the reference position
        # No action is taken
      }
    }
    return(list(base = c_base, qual = c_qual, indel = c_indel))
  }
}


#' Get the informations about the presence of a insertion in the read.
#'
#' @inheritParams build_fragments_info_table*
#' @param pos numeric value representing the position of interest
#' @param alt character vector representing the insertion sequence
#' @param r_pos numeric value representing the read mapping position
#' @param r_cigar character vector representing the read CIGAR
#' @param r_query character vector representing the read base sequence
#' @param r_qual character vector representing the read sequencing qualities
#' @return a named list with names `base`,`qual` and `indel`.
#'
#' @importFrom stringr str_extract utils
#'
#' @noRd
get_insertion <- function(pos, alt, r_pos, r_cigar, r_query, r_qual) {
  c_base <- NA
  c_qual <- NA
  c_indel <- "No_support_ins"  # By default, no insertion
  alt_len <- nchar(alt)

  # Parse the CIGAR string
  ops <- parse_cigar(r_cigar)

  # Iterating over CIGAR operations
  current_pos <- r_pos  # Current reference position, aligned with the read
  for (i in seq_len(nrow(ops))) {
    op_len <- ops$length[i]
    op_type <- ops$type[i]

    if (op_type == "M") {
      # M: match or mismatch (but advances the reference position)
      current_pos <- current_pos + op_len
    } else if (op_type == "D") {
      # D: deletion in the read (advances the reference position but not the read)
      # Therefore, we modify current_pos
      current_pos <- current_pos + op_len
    } else if (op_type == "I") {
      # I: insertion in the read (addition in the read compared to the reference)
      # Compare this insertion with the reference
      diff_mutation_bam <- pos - (current_pos - 1)
      is_length_ok <- (alt_len == op_len)

      ins_rep <- define_ins_rep(alt_len, current_pos, r_pos, r_cigar, r_query)

      c_indel <- check_seq_rep(ins_rep, is_length_ok, diff_mutation_bam, op_len)

      if (c_indel == 1) {
        c_qual <- find_c_qual(alt_len, current_pos, r_pos, r_cigar, r_qual)
        # A corresponding insertion was found
        break
      }
    } else if (op_type %in% c("N", "=", "X")) {
      # N: skip in the reference, =: perfect match, X: mismatch
      # All advance the reference position
      current_pos <- current_pos + op_len
    } else {
      # S, H, P: Soft clip, Hard clip, Pad - do not advance the reference position
      # No action is taken
    }
  }

  return(list(base = c_base, qual = c_qual, indel = c_indel))
}





#' Get genomic positions of deletions and insertions of reads.
#'
#' @inheritParams get_nb_unique_del_ins
#' @param r_pos numeric value representing the read mapping position
#' @param r_cigar character vector representing the read CIGAR
#' @return a named list with names `d_pos` and `i_pos`.
#'
#' @importFrom stringr str_extract
#'
#' @noRd
get_pos_indel_from_read <- function(r_pos, r_cigar){
  # table of operations-cursor movement
  Op <- c("M", "I", "D", "N", "S", "H", "P", "=", "X")
  Consumes_query <- c("yes", "yes", "no", "no", "yes", "no", "no", "yes", "yes")
  Consumes_reference <- c("yes", "no", "yes", "yes", "no", "no", "no", "yes", "yes")
  CIGAR <- data.frame(Op, Consumes_query, Consumes_reference)


  # cursors and index to be incremented
  cursor_reference <- r_pos
  d_pos <- c()
  i_pos <- c()

  while (nchar(r_cigar)>0){
    # get current cigar operation
    c_dg_str <- str_extract(r_cigar, "^[\\d]+")
    c_dg_num <- as.numeric(c_dg_str)
    c_op <- str_extract(r_cigar, "(?<=[\\d])[MIDNSHP=X]")
    r_cigar_left <- gsub(paste0("^", c_dg_str, c_op), "", r_cigar)


    # execute cigar operation to move cursors
    c_reference <- CIGAR[CIGAR$Op==c_op, "Consumes_reference"]


    # register position if deletion or insertion
    if (c_op=="D"){
      d_pos <- append(d_pos, cursor_reference)
    } else if (c_op=="I"){
      i_pos <- append(i_pos, cursor_reference )
    }

    # execute r_cigar operation along the reference
    if (c_reference=="yes"){
      cursor_reference <- cursor_reference + 1
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

  list(d_pos=d_pos, i_pos=i_pos)
}


#' Get the count of unique deletions and insertions.
#'
#' @param r_pos1 numeric value representing the read 1 mapping position
#' @param r_pos2 numeric value representing the read 2 mapping position
#' @param r_cigar1 character vector representing the read 1 CIGAR
#' @param r_cigar2 character vector representing the read 2 CIGAR
#' @return a named list with names `n_unique_del` and `n_unique_ins`.
#'
#' @importFrom stringr str_extract
#'
#' @noRd
get_nb_unique_del_ins <- function(r_pos1, r_cigar1, r_pos2, r_cigar2){

     # genomic positions of deletions and insertions from read 1
     d_pos_i_pos1 <- get_pos_indel_from_read(r_pos=r_pos1, r_cigar=r_cigar1)
     d_pos1 <- d_pos_i_pos1$d_pos
     i_pos1 <- d_pos_i_pos1$i_pos

     # genomic positions of deletions and insertions from read 2
     d_pos_i_pos2 <- get_pos_indel_from_read(r_pos=r_pos2, r_cigar=r_cigar2)
     d_pos2 <- d_pos_i_pos2$d_pos
     i_pos2 <- d_pos_i_pos2$i_pos

     # counts of unique deletions
     d_pos1_setd <- setdiff(d_pos1, d_pos2) #value that exists only in the first vector d_pos1
     d_pos2_setd <- setdiff(d_pos2, d_pos1) #value that exists only in the first vector d_pos2
     d_pos_int <- intersect(d_pos1, d_pos2) #identical values between the two vectors d_pos1 and d_pos2
     d_pos1_id <- identical(d_pos_int, d_pos1) #if the two vectors d_pos_int and d_pos1 are equal then this returns TRUE otherwise it returns FALSE
     d_pos2_id <- identical(d_pos_int, d_pos2) #if the two vectors d_pos_int and d_pos2 are equal then this returns TRUE otherwise it returns FALSE
     if (d_pos1_id == TRUE | d_pos2_id == TRUE){
       n_unique_del <- length(d_pos1_setd) + length(d_pos2_setd) + length(d_pos_int)
     } else {
        # print(d_pos1)
        # print(d_pos2)
       n_unique_del <- length(d_pos1) + length(d_pos2)
     }

     # counts of unique insertions
     i_pos_12 <- list(i_pos1, i_pos2)
     i_pos_int <- intersect(i_pos1[1], i_pos2[1]) #identical values between the first two vectors i_pos1 and i_pos2
     i_pos1_id <- identical(i_pos_int, i_pos1[1]) #if the two vectors i_pos_int and i_pos1 are equal then this returns TRUE otherwise it returns FALSE
     i_pos2_id <- identical(i_pos_int, i_pos2[1]) #if the two vectors i_pos_int and i_pos2 are equal then this returns TRUE otherwise it returns FALSE
     if (i_pos1_id == TRUE & i_pos2_id == TRUE){
       n_unique_ins <- max(lengths(i_pos_12))
     } else {
        # print(i_pos1)
        # print(i_pos2)
       n_unique_ins <- length(i_pos1) + length(i_pos2)
     }

  list(unique_del=n_unique_del, unique_ins=n_unique_ins)
}

#' Get inner distance between two reads from the same fragment.
#'
#' @param r_pos1 numeric value representing the read 1 mapping position
#' @param r_pos2 numeric value representing the read 2 mapping position
#' @param r_cigar1 character vector representing the read 1 CIGAR
#' @param r_cigar2 character vector representing the read 2 CIGAR
#' @param read_length1 numeric value representing the length of read 1
#' @param read_length2 numeric value representing the length of read 2
#' @return a named list with names `a` and `b`.
#'
#' @importFrom stringr str_extract
#'
#' @noRd
get_fragment_inner_distance <- function(r_pos1, r_cigar1, read_length1, r_pos2, r_cigar2, read_length2){
    # bases matched read 1/2
    bases_match1 <- sum(as.numeric(str_extract_all(r_cigar1, "[[:digit:]]+(?=M)", simplify=TRUE)))
    bases_match2 <- sum(as.numeric(str_extract_all(r_cigar2, "[[:digit:]]+(?=M)", simplify=TRUE)))

    # bases deleted read 1/2
    bases_del1 <- sum(as.numeric(str_extract_all(r_cigar1, "[[:digit:]]+(?=D)", simplify=TRUE)))
    bases_del2 <- sum(as.numeric(str_extract_all(r_cigar2, "[[:digit:]]+(?=D)", simplify=TRUE)))

    # bases inserted read 1/2
    bases_ins1 <- sum(as.numeric(str_extract_all(r_cigar1, "[[:digit:]]+(?=I)", simplify=TRUE)))
    bases_ins2 <- sum(as.numeric(str_extract_all(r_cigar2, "[[:digit:]]+(?=I)", simplify=TRUE)))

    # bases soft-clipped left read 1/2
    bases_softcl_left1 <- as.numeric(str_extract(r_cigar1, "^[[:digit:]]+(?=S)"))
    bases_softcl_left1 <- ifelse(is.na(bases_softcl_left1), 0, bases_softcl_left1)
    bases_softcl_left2 <- as.numeric(str_extract(r_cigar2, "^[[:digit:]]+(?=S)"))
    bases_softcl_left2 <- ifelse(is.na(bases_softcl_left2), 0, bases_softcl_left2)

    # bases soft-clipped right read 1/2
    bases_softcl_right1 <- as.numeric(str_extract(r_cigar1, "[[:digit:]]+(?=S$)"))
    bases_softcl_right1 <- ifelse(is.na(bases_softcl_right1), 0, bases_softcl_right1)
    bases_softcl_right2 <- as.numeric(str_extract(r_cigar2, "[[:digit:]]+(?=S$)"))
    bases_softcl_right2 <- ifelse(is.na(bases_softcl_right2), 0, bases_softcl_right2)

    # number of soft clipped bases genomically left
    if (r_pos1 <= r_pos2){
      bases_softcl_left <- bases_softcl_left1
      bases_softcl_right <- bases_softcl_right2
    } else {
      bases_softcl_left <- bases_softcl_left2
      bases_softcl_right <- bases_softcl_right1
    }


    # inner distance (two formulas)
    if (r_pos1 <= r_pos2){
      # 1 = LEFT
      # 2 = RIGHT
      inner_distance_a <- r_pos2 - r_pos1 + bases_softcl_left1 - read_length1 - bases_softcl_left2 + bases_ins1 - bases_del1

      inner_distance_b <- r_pos2 - r_pos1 - bases_match1 - bases_del1 - bases_softcl_left2 - bases_softcl_right1
    } else {
      # 2 = LEFT
      # 1 = RIGHT
      inner_distance_a <- r_pos1 - r_pos2 + bases_softcl_left2 - read_length2 - bases_softcl_left1 + bases_ins2 - bases_del2

      inner_distance_b <- r_pos1 - r_pos2 - bases_match2 - bases_del2 - bases_softcl_left1 - bases_softcl_right2
    }

  list(a=inner_distance_a, b=inner_distance_b)
}



#' Prepare R sockets for running in parallel.
#'
#' @inheritParams build_fragments_info_table
#' @return an object of class "cluster"
#'
#' @importFrom parallel makePSOCKcluster clusterExport clusterEvalQ
#' @importFrom doParallel registerDoParallel
#'
#' @noRd
setup_parallel_computations <- function(n_cores){
  cl <- makePSOCKcluster(n_cores)
  registerDoParallel(cl)
  # tmp <- clusterEvalQ(cl, library(fRagmentomics))
  cl
}


#' Create a table with statistics about fragments covering a position of interest (e.g the position of a mutation).
#'
#' @param df_sam dataframe containing reads from SAM format
#' @param chr chromosome for position of interest
#' @param pos chromosome for position of interest
#' @param ref reference sequence for position of interest
#' @param alt alternative sequence for position of interest
#' @param del_info character vector in the format "$pos_final,$rep_del"
#' @param mutation_type type of the mutation in insertion, deletion and mutation
#' @param n_cores the number of cores available
#' @return a dataframe
#'
#' @importFrom readr read_tsv
#' @importFrom stringr str_extract_all
#' @importFrom parallel stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom utils txtProgressBar setTxtProgressBar flush.console
#' @export
build_fragments_info_table <- function(df_sam, sample_id, chr, pos, ref, alt, del_info, mutation_type, n_cores=1){

  # list of all fragments names
  fragments_names <- unique(df_sam[,1,drop=T])
  n_fragments <- length(fragments_names)

  # prepare cluster
  cat(paste("-exporting functions to cluster..."))
  cl <- setup_parallel_computations(n_cores)
  cat("done!\n")

  # function to rbind with progressbar
  rbindprogress <- function(iterator){
    pb <- txtProgressBar(min = 1, max = iterator - 1, style = 3)
    count <- 0
    function(...) {
      count <<- count + length(list(...)) - 1
      setTxtProgressBar(pb, count)
      flush.console()
      rbind(...)
    }
  }

  # run a loop on fragments names
  i <- NULL
  df_fragments_info <- foreach(i=1:n_fragments,
                               .combine = rbindprogress(n_fragments),
                               .inorder = F,
                               .multicombine = F) %dopar% {
    fragment_name <- fragments_names[[i]]

    # select reads corresponding to fragment
    df_fragment_reads <- df_sam[df_sam[,1,drop=T]==fragment_name,,drop=F]

    # sanity checks
    fragment_qc <- ""

    if (!all(df_fragment_reads[,3]==chr)){
      fragment_qc <- paste(fragment_qc, "&", paste("read(s) from other chromosome than", chr))
    }

    if (nrow(df_fragment_reads)!=2){
      fragment_qc <- paste(fragment_qc, "&", paste(nrow(df_fragment_reads), "read(s)"))
    }

    if (fragment_qc != ""){
      # fragments does not pass checks
      df_fragment <- data.frame(Chromosome=chr, Position=pos, Fragment=fragment_name, Fragment_Check=fragment_qc)
    } else {
      ####
      #### load reads 1 and 2 into one-line dataframes
      ####

      # read 1 is defined as the read having 0x40 in its flag while read 2 has 0x80 in its flag value
      flag_a <- df_fragment_reads[1,"flag",drop=T]
      flag_b <- df_fragment_reads[2,"flag",drop=T]

      flag_a <- bitvalues_from_bam_flag(as.integer(flag_a), bitnames=c("isFirstMateRead","isSecondMateRead"))
      flag_b <- bitvalues_from_bam_flag(as.integer(flag_b), bitnames=c("isFirstMateRead","isSecondMateRead"))

      if ( (flag_a[,"isFirstMateRead"] == 1) & (flag_b[,"isSecondMateRead"] == 1) ){
        df_read1 <- df_fragment_reads[1, ]
        df_read2 <- df_fragment_reads[2, ]
      } else if ( (flag_b[,"isFirstMateRead"] == 1) & (flag_a[,"isSecondMateRead"] == 1) ){
        df_read1 <- df_fragment_reads[2, ]
        df_read2 <- df_fragment_reads[1, ]
      } else {
        stop("invalid read flags. One read should be first mate, the other second mate:\n", flag_a, "\n", flag_b)
      }

      ####
      #### retrieve read-level stats
      ####

      # mapping quality
      mapq1 <- df_read1$mapq
      mapq2 <- df_read2$mapq

      # tlen
      tlen1 <- df_read1$tlen
      tlen2 <- df_read2$tlen

      # reads cigar
      r_cigar1 <- df_read1$cigar
      r_cigar2 <- df_read2$cigar

      # reads mapped positions
      r_pos1 <- df_read1$pos
      r_pos2 <- df_read2$pos

      # reads queries
      r_query1 <- df_read1$seq
      r_query2 <- df_read2$seq

      # reads qualities
      r_qual1 <- df_read1$qual
      r_qual2 <- df_read2$qual

      # reads length
      read_length1 <- nchar(r_query1)
      read_length2 <- nchar(r_query2)

      # identity variable to identify the position of read 1 and 2 depending on whether they are on the left or right
      if (r_pos1 <= r_pos2){
        identity_RIGHT <- 2
        identity_LEFT <- 1
      } else {
        identity_RIGHT <- 1
        identity_LEFT <- 2
      }

      ####
      #### retrieve base and qualities
      ####

      if (mutation_type == "mutation") {
        r_info1 <- get_base_qual_from_read(pos=pos, r_pos=r_pos1, r_cigar=r_cigar1, r_query=r_query1, r_qual=r_qual1)
        r_info2 <- get_base_qual_from_read(pos=pos, r_pos=r_pos2, r_cigar=r_cigar2, r_query=r_query2, r_qual=r_qual2)

      } else if (mutation_type == "deletion") {
        r_info1 <- get_deletion(pos=pos, ref=ref, del_info=del_info, r_pos=r_pos1, r_cigar=r_cigar1)
        r_info2 <- get_deletion(pos=pos, ref=ref, del_info=del_info, r_pos=r_pos2, r_cigar=r_cigar2)

      } else if (mutation_type == "insertion") {
        r_info1 <- get_insertion(pos=pos, alt=alt, r_pos=r_pos1, r_cigar=r_cigar1, r_query=r_query1, r_qual=r_qual1)
        r_info2 <- get_insertion(pos=pos, alt=alt, r_pos=r_pos2, r_cigar=r_cigar2, r_query=r_query2, r_qual=r_qual2)

      } else {
        r_info1 <- "Mutation_type is not defined. You should define it as 'mutation', 'deletion' or 'insertion'"
        r_info2 <- "Mutation_type is not defined. You should define it as 'mutation', 'deletion' or 'insertion'"
      }


      ####
      #### retrieve bases : inserted, deleted, soft-clipped left and soft-clipped right
      ####

      # bases soft-clipped left read 1/2
      bases_softcl_left1 <- as.numeric(str_extract(r_cigar1, "^[[:digit:]]+(?=S)"))
      bases_softcl_left1 <- ifelse(is.na(bases_softcl_left1), 0, bases_softcl_left1)
      bases_softcl_left2 <- as.numeric(str_extract(r_cigar2, "^[[:digit:]]+(?=S)"))
      bases_softcl_left2 <- ifelse(is.na(bases_softcl_left2), 0, bases_softcl_left2)

      # bases soft-clipped right read 1/2
      bases_softcl_right1 <- as.numeric(str_extract(r_cigar1, "[[:digit:]]+(?=S$)"))
      bases_softcl_right1 <- ifelse(is.na(bases_softcl_right1), 0, bases_softcl_right1)
      bases_softcl_right2 <- as.numeric(str_extract(r_cigar2, "[[:digit:]]+(?=S$)"))
      bases_softcl_right2 <- ifelse(is.na(bases_softcl_right2), 0, bases_softcl_right2)

      # number of soft clipped bases genomically left
      if (r_pos1 <= r_pos2){
        bases_softcl_left <- bases_softcl_left1
        bases_softcl_right <- bases_softcl_right2
      } else {
        bases_softcl_left <- bases_softcl_left2
        bases_softcl_right <- bases_softcl_right1
      }

      # counts of unique deletions and insertions
      n_unique_del_ins <- get_nb_unique_del_ins(r_pos1, r_cigar1, r_pos2, r_cigar2)
      n_unique_del <- n_unique_del_ins$unique_del
      n_unique_ins <- n_unique_del_ins$unique_ins

      ####
      #### retrieve inner distance and absolute_size
      ####
      inner_distances <- get_fragment_inner_distance(r_pos1, r_cigar1, read_length1, r_pos2, r_cigar2, read_length2)
      inner_distance_a <- inner_distances$a
      inner_distance_b <- inner_distances$b

      # # Check that the two formulas for calculating inner distance return the same values
      # if (inner_distance_a == inner_distance_b){
      #   inner_distance <- inner_distance_a
      # } else {
      #   inner_distance <- "error"
      # }

      # absolute size
      absolute_size <- read_length1 + read_length2 + inner_distance_b

      # result dataframe
      df_fragment <- data.frame(Sample_Id=sample_id, Chromosome=chr, Position=pos, Ref=ref, Alt=alt, Fragment=fragment_name, Fragment_Check="OK",
                                Mutation_type=mutation_type, Mapq_Read_1=mapq1, Mapq_Read_2=mapq2, Base_Read_1=r_info1$base,
                                Base_Read_2=r_info2$base, Qual_Read_1=r_info1$qual, Qual_Read_2=r_info2$qual,
                                Indel_1=r_info1$indel, Indel_2=r_info2$indel,
                                TLEN=abs(tlen1), Read_length_1=read_length1, Read_length_2=read_length2,
                                Startpos_read_1=r_pos1, Startpos_read_2=r_pos2,
                                Identity_RIGHT=identity_RIGHT, Identity_LEFT=identity_LEFT,
                                Bases_del=n_unique_del, Bases_ins=n_unique_ins,
                                Bases_soft_clip_left=bases_softcl_left, Bases_soft_clip_right=bases_softcl_right,
                                Inner_distance_a=inner_distance_a, Inner_distance_b=inner_distance_b,
                                Absolute_size=absolute_size)
    }
    df_fragment
  }

  # terminate cluster
  stopCluster(cl)

  df_fragments_info
}
