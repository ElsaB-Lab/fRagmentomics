#' Retrieve a genomic sequence from a FASTA source
#'
#' @description Fetches a DNA sequence for a specified genomic region. This function can
#' operate in two modes: either by reading directly from an indexed FASTA file
#' or by extracting a subsequence from a pre-loaded sequence object.
#'
#' @inheritParams normalize_to_vcf_rep
#' @param start Start position to be retrieved
#' @param end End position to be retrieved
#' @param fasta_seq A list with the fasta sequence between two positions.
#'
#' @return A character string representing the single nucleotide at
#' the specified chromosome and position.
#'
#' @importFrom Biostrings getSeq
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#'
#' @noRd
get_seq_from_fasta <- function(chr, start, end, fasta_fafile = NULL, fasta_seq = NULL) {
    if (!is.null(fasta_seq)) {
        f_chr <- fasta_seq$chr
        f_start <- fasta_seq$start
        f_end <- fasta_seq$end
        f_seq <- fasta_seq$seq
        if (f_chr != chr) {
            stop(sprintf(
                "Requested chromosome '%s' does not match the available chromosome '%s'.",
                chr,
                f_chr
            ))
        }
        if (start < f_start || end > f_end) {
            stop(sprintf(
                "Requested sequence range %d:%d does not fit into the available reference sequence range %d:%d.",
                start,
                end,
                f_start,
                f_end
            ))
        }
        ref_seq <- substr(f_seq, start - f_start + 1, end - f_start + 1)
    } else {
        # Fetch a single nucleotide from the FASTA
        # using the chromosome, and start = end = position
        ref_seq <- Biostrings::getSeq(
            x = fasta_fafile,
            param = GenomicRanges::GRanges(
                seqnames = chr,
                ranges = IRanges::IRanges(start = start, end = end)
            )
        )
    }

    # Transform fasta_seq into string
    ref_seq_char <- unname(as.character(ref_seq))
    return(as.character(ref_seq_char))
}


#' Harmonize chromosome notation to match a FASTA file
#'
#' @description This utility function adjusts a chromosome name to match the style used in a reference FASTA file's
#' index (e.g., converting "1" to "chr1" or vice versa).
#'
#' @inheritParams normalize_to_vcf_rep
#' @importFrom Rsamtools scanFaIndex
#' @importFrom GenomeInfoDb seqnames
#'
#' @return String with the chr notation harmonized to match the FASTA reference.
#'
#' @noRd
harmonize_chr_to_fasta <- function(chr, fasta_fafile) {
    # Extract chromosome names from the FASTA index
    fasta_index <- scanFaIndex(fasta_fafile) # Get the indexed FASTA sequence info
    fasta_chromosomes <- seqnames(fasta_index) # Extract chromosome names

    # Determine whether FASTA uses "chr1" or "1"
    fasta_format <- ifelse(any(grepl("^chr", fasta_chromosomes)), "chr", "no_chr")

    # Harmonize chromosome notation
    if (fasta_format == "chr" && !grepl("^chr", chr)) {
        chr <- paste0("chr", chr) # Add "chr" if missing
    } else if (fasta_format == "no_chr" && grepl("^chr", chr)) {
        chr <- sub("^chr", "", chr) # Remove "chr" if present
    }

    chr
}


#' Check that the reference allele matches with the FASTA sequence
#'
#' @description This function checks whether the reference allele ('ref') at a given genomic
#' position ('pos') matches the corresponding nucleotide in the provided FASTA
#' reference genome.
#'
#' @inheritParams normalize_to_vcf_rep
#'
#' @return Logical 'TRUE' if the reference allele matches the FASTA nucleotide,
#'   otherwise 'FALSE' with a warning.
#'
#' @importFrom Biostrings getSeq
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
#' @noRd
check_if_ref_matches_fasta <- function(chr, pos, ref, fasta_fafile) {
    # Fetch the names (i.e., FASTA headers) for each sequence
    idx <- scanFaIndex(fasta_fafile)
    seq_names <- as.character(seqnames(idx))

    # Check if the chrom existe in the Fasta
    if (!chr %in% seq_names) {
        warning(sprintf("Chromosome '%s' not found in FASTA.", chr))
        return(FALSE)
    }

    # Calculate the end position in the FASTA based on the length of 'ref'
    ref_length <- nchar(ref)
    end_pos <- pos + ref_length - 1

    # Retrieve the expected reference sequence from FASTA
    fasta_seq <- get_seq_from_fasta(chr = chr, start = pos, end = end_pos, fasta_fafile = fasta_fafile)

    # Compare the full 'ref' to the fetched FASTA sequence
    if (identical(ref, fasta_seq)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}
