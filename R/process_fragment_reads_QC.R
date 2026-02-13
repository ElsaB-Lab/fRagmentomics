#' Perform quality control (QC) checks on the reads of a single fragment
#'
#' @description This function runs a series of validation checks on the reads belonging to a single DNA fragment.
#' It ensures the pair is complete, properly mapped, and topologically consistent.
#'
#' @details
#' The logic handles the following scenarios:
#' \itemize{
#'   \item **1 Read**: If a read is missing, it deduces why (Unmapped, Translocation, or Far Away).
#'   \item **2 Reads**: It verifies:
#'     \itemize{
#'       \item *Mapping*: Both reads must have valid positions (POS != 0).
#'       \item *Chromosome*: Both reads must be on the expected chromosome.
#'       \item *Consistency*: The mate information (RNEXT) must match the current chromosome.
#'     }
#' }
#'
#' @inheritParams extract_fragment_features
#' @param df_fragment_reads a data frame containing reads for a single fragment.
#'
#' @return a single string containing all quality control messages. Returns "" if OK.
#'
#' @noRd
process_fragment_reads_qc <- function(df_fragment_reads, chr) {
    qc_messages <- c()
    num_reads <- nrow(df_fragment_reads)

    # --- CASE A: Abnormal Read Count (Not 2) ---
    if (num_reads != 2) {
        # Detailed diagnosis if we have exactly 1 read
        if (num_reads == 1) {
            msg <- paste("Fragment has 1 read")
            read <- df_fragment_reads[1, , drop = FALSE]

            # Normalize mate chromosome info
            mate_chr <- read$RNEXT
            if (!is.na(mate_chr) && mate_chr == "=") {
                mate_chr <- read$RNAME
            }

            # Diagnosis 1: Mate is unmapped (RNEXT is *)
            if (is.na(mate_chr) || mate_chr == "*") {
                msg <- paste(msg, "- Mate is unmapped")
            }
            # Diagnosis 2: Translocation (Mate on diff chr, filtered out by read_bam)
            else if (mate_chr != read$RNAME) {
                msg <- paste(msg, "- Mate maps to a different chromosome:", mate_chr)
            }
            # Diagnosis 3: Far Away (Mate on same chr but outside window, filtered out)
            else if (mate_chr == read$RNAME) {
                dist <- if (!is.na(read$TLEN)) abs(read$TLEN) else "NA"
                msg <- paste(msg, "- Mate maps outside loaded region. TLEN:", dist)
            }
        } else {
            msg <- paste("Fragment has", num_reads, "read(s)")
        }

        qc_messages <- c(qc_messages, msg)
    }

    # --- CASE B: Normal Read Count (2 Reads) - Check Consistency ---
    if (num_reads == 2) {
        # Use first read to check pairing info consistency
        read <- df_fragment_reads[1, , drop = FALSE]

        # Normalize mate chromosome
        mate_chr <- read$RNEXT
        if (!is.na(mate_chr) && mate_chr == "=") {
            mate_chr <- read$RNAME
        }

        # Test 1: Wrong Chromosome (Are the loaded reads actually on the target chr?)
        if (!all(df_fragment_reads$RNAME == chr)) {
            qc_messages <- c(qc_messages, paste("Read(s) found on a chromosome other than", chr))
        }

        # Test 2: Missing Mate Info (RNEXT is *)
        if (is.na(read$RNEXT) || read$RNEXT == "*") {
            qc_messages <- c(qc_messages, "Mate chromosome info is not available")
        }

        # Test 3: Internal Consistency (Do reads claim to be on different chromosomes?)
        # Even if both reads are present, RNEXT should match RNAME.
        else if (mate_chr != read$RNAME) {
            qc_messages <- c(qc_messages, paste("- Mate maps to a different chromosome:", mate_chr))
        }

        # Test 4: Mapping Status
        # Verify that POS and PNEXT are strictly positive (mapped)
        if (any(is.na(df_fragment_reads$POS)) || any(df_fragment_reads$POS == 0)) {
            qc_messages <- c(qc_messages, "One or both reads are unmapped (POS=0 or NA)")
        }
    }

    fragment_qc <- paste(qc_messages, collapse = " & ")
    return(fragment_qc)
}
