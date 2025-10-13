#' Compute fragment size
#'
#' @description
#' Calculates the size of a DNA fragment from its aligned paired-end reads. The
#' function provides a more precise estimate than the SAM/BAM TLEN field by
#' accounting for complex alignment features like soft-clipping and indels not
#' only based one genomic position.
#'
#' @param read_stats_5p 5p read infos
#' @param read_stats_3p 3p read infos
#'
#' @return a integer corresponding to the fragment size
#'
#' @importFrom stringr str_extract str_extract_all
#'
#' @keywords internal
get_fragment_size <- function(read_stats_5p, read_stats_3p) {
    # --- Extract necessary metrics --- bases matched read
    bases_match_5p <- sum(as.numeric(
        stringr::str_extract_all(read_stats_5p$CIGAR,
        "[[:digit:]]+(?=M)", simplify = TRUE)))

    # bases deleted read
    bases_del_5p <- sum(as.numeric(
        stringr::str_extract_all(read_stats_5p$CIGAR,
        "[[:digit:]]+(?=D)", simplify = TRUE)))

    # bases soft-clipped left
    bases_softcl_left_3p <- as.numeric(
        str_extract(read_stats_3p$CIGAR, "^[[:digit:]]+(?=S)"))
    bases_softcl_left_3p <- ifelse(
        is.na(bases_softcl_left_3p), 0, bases_softcl_left_3p)

    # bases soft-clipped right
    bases_softcl_right_5p <- as.numeric(
        str_extract(read_stats_5p$CIGAR, "[[:digit:]]+(?=S$)"))
    bases_softcl_right_5p <- ifelse(
        is.na(bases_softcl_right_5p), 0, bases_softcl_right_5p)

    # --- Define overlapping windows ---
    start_overlap <- read_stats_3p$POS - bases_softcl_left_3p
    end_overlap <- read_stats_5p$POS + bases_match_5p + bases_del_5p +
        bases_softcl_right_5p - 1

    # --- Calculate inner size of the fragment ---
    inner_distance <- start_overlap - end_overlap - 1

    # Taking indels into account in the overlapping windows to no count them
    # twice Get indels for each read
    indels_5p <- get_pos_indels_from_read(read_stats_5p)
    indels_3p <- get_pos_indels_from_read(read_stats_3p)

    # Merge the deletion lists and keep unique positions in the overlapping
    # windows
    unique_fragment_deletions <- unique(
        c(indels_5p$deletions, indels_3p$deletions))
    mask_in_overlap <- unique_fragment_deletions >= start_overlap &
        unique_fragment_deletions <= end_overlap
    overlapping_window_deletion <- unique_fragment_deletions[mask_in_overlap]

    # Merge the insertion lists and keep unique positions in the overlapping
    # windows
    unique_fragment_insertions <- unique(
        c(indels_5p$insertions, indels_3p$insertions))
    overlapping_window_insertion <- numeric(0)
    if (length(unique_fragment_insertions) > 0) {
        mask_in_overlap <- unique_fragment_insertions >= start_overlap &
            unique_fragment_insertions <= end_overlap
        overlapping_window_insertion <-
            unique_fragment_insertions[mask_in_overlap]
    }

    fragment_size <- read_stats_5p$read_length + read_stats_3p$read_length +
        inner_distance + length(overlapping_window_deletion) -
        length(overlapping_window_insertion)
}
