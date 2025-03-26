# Project : ElsaBLab_fRagmentomics

#' Sanity check: Validate reference allele against FASTA sequence
#'
#' This function checks whether the reference allele (`ref`) at a given genomic
#' position (`pos`) matches the corresponding nucleotide in the provided FASTA
#' reference genome.
#'
#' @inheritParams process_fragment
#' @param fasta A FASTA ref genome (loaded using `Biostrings` or a file path).
#'
#' @return Logical `TRUE` if the reference allele matches the FASTA nucleotide,
#'   otherwise `FALSE` with a warning.
#'
#' @importFrom Biostrings getSeq
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
#' @noRd
sanity_check_ref_in_fasta <- function(chr, pos, ref, fasta) {
  # Fetch the names (i.e., FASTA headers) for each sequence
  idx <- scanFaIndex(fasta)
  seq_names <- as.character(seqnames(idx))

  # Check if the chrom existe in the Fasta
  if (!chr %in% seq_names) {
    warning(paste0("Chromosome ", chr, " not found in FASTA."))
    return(FALSE)
  }

  # Calculate the end position in the FASTA based on the length of `ref`
  ref_length <- nchar(ref)
  end_pos <- pos + ref_length - 1

  # Retrieve the expected reference sequence from FASTA
  fasta_seq <- Biostrings::getSeq(
    x = fasta,
    param = GRanges(
      seqnames = chr,
      ranges = IRanges(start = pos, end = end_pos)
    )
  )

  # Transform fasta_seq into string
  fasta_seq_char <- unname(as.character(fasta_seq))

  # Compare the full `ref` to the fetched FASTA sequence
  if (identical(ref, fasta_seq_char)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
