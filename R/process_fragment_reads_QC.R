#' Perform quality control (QC) checks on the reads of a single fragment
#'
#' @description This function runs a series of validation checks on the reads belonging to a single DNA fragment to
#' ensure they form a valid, properly-mapped pair.
#'
#' @details
#' The function performs the following specific QC checks:
#' \itemize{
#'   \item **Paired Read Count**: Verifies that the fragment consists of exactly two reads.
#'   \item **Chromosome Consistency**: Ensures all reads in the fragment map to the expected chromosome ('chr').
#'   \item **Mate Information**: Checks that the mate-pair chromosome ('RNEXT') is available (not '*').
#'   \item **Intra-chromosomal Mapping**: Confirms that the mate read maps to the same chromosome.
#'   \item **Mapping Status**: Verifies that both the read and its mate are mapped (i.e., 'POS' and 'PNEXT' are not 0 or 'NA').
#' }
#'
#' @inheritParams extract_fragment_features
#' @param df_fragment_reads a data frame containing all reads for a single fragment.
#'
#' @return a single string containing all quality control messages.
#'
#' @noRd
process_fragment_reads_qc <- function(df_fragment_reads, chr) {
    qc_messages <- c()

    # Test 1: Report if the fragment does not consist of exactly two reads.
    if (nrow(df_fragment_reads) != 2) {
        qc_messages <- c(qc_messages, paste("Fragment has", nrow(df_fragment_reads), "read(s)"))
    }

    # Test 2: Report if any read in the fragment is not on the expected chromosome.
    if (!all(df_fragment_reads$RNAME == chr)) {
        qc_messages <- c(qc_messages, paste("Read(s) found on a chromosome other than", chr))
    }

    if (nrow(df_fragment_reads) > 0) {
        # Isolate the first read of the fragment
        read <- df_fragment_reads[1, , drop = FALSE]

        # Test 3: Report if mate chromosome info is missing (* in BAM).
        if (is.na(read$RNEXT)) {
            qc_messages <- c(qc_messages, "Mate chromosome info is not available (RNEXT was *)")
        }

        # Test 4: Report if the mate appears to be on another chromosome (translocation).
        else if (read$RNEXT != "=" && read$RNAME != read$RNEXT) {
            qc_messages <- c(qc_messages, paste("Mate maps to a different chromosome. RNEXT=", read$RNEXT))
        }

        # Test 5: Report if the read or its mate is marked as unmapped.
        if (is.na(read$POS) || read$POS == 0 || is.na(read$PNEXT) || read$PNEXT == 0) {
            qc_messages <- c(qc_messages, paste("Read or mate is unmapped. POS=", read$POS, "& PNEXT=", read$PNEXT))
        }
    }
    fragment_qc <- paste(qc_messages, collapse = " & ")

    return(fragment_qc)
}
