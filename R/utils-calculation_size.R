#' @title Compute fragment inner distance
#' @description Calculates the inner distance between two paired-end sequencing reads.
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