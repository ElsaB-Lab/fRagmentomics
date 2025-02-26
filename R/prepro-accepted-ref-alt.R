#' Retrieve a single nucleotide from a multi-chromosome FASTA
#'
#' This function retrieves a single nucleotide from a multi-chromosome FASTA reference
#' at the specified chromosome and position (1-based coordinates).
#'
#' @param chr A character string specifying the chromosome name "chr1".
#' @param position An integer specifying the genomic position (1-based).
#' @param fasta A FASTA reference that contains all chromosomes. Any reference object
#'   compatible with Biostrings::getSeq.
#'
#' @return A character string representing the single nucleotide at the specified
#'   chromosome and position.
#' 
#' @importFrom Biostrings getSeq
retrieve_anchor_base <- function(chr, position, fasta) {
  # Fetch a single nucleotide from the FASTA
  # using the chromosome, and start = end = position
  anchor_base <- as.character(
    Biostrings::getSeq(
      x     = fasta,
      names = chr,
      start = position,
      end   = position
    )
  )
  return(anchor_base)
}

#' Trim Common Left Prefix from Two Sequences
#'
#' This function trims any common left prefix between two sequences while ensuring
#' that at least one shared base remains. If a common prefix is found, the function
#' keeps exactly one shared base and removes the rest of the overlapping prefix.
#'
#' @param a A character string representing the first sequence.
#' @param b A character string representing the second sequence.
#'
#' @return A list with two elements:
#'   The modified first sequence with minimal shared prefix.}
#'   The modified second sequence with minimal shared prefix.}
#'
#' @importFrom stringr str_sub
trim_left_prefix <- function(a, b) {
  min_len <- min(nchar(a), nchar(b))
  prefix_length <- 0
  
  while (prefix_length < min_len &&
         substr(a, prefix_length + 1, prefix_length + 1) ==
         substr(b, prefix_length + 1, prefix_length + 1)) {
    prefix_length <- prefix_length + 1
  }
  
  # Keep exactly one shared base if a prefix exists
  if (prefix_length > 0) {
    prefix_length <- prefix_length - 1
  }
  
  new_a <- substr(a, prefix_length + 1, nchar(a))
  new_b <- substr(b, prefix_length + 1, nchar(b))
  
  return(list(a = new_a, b = new_b))
}

#' Normalize REF and ALT to standard VCF conventions
#'
#' This function normalizes REF and ALT alleles to the standard VCF conventions.
#' It handles both the common VCF-like formats (SNV, MNV, insertion, deletion)
#' and also the alternate insertion/deletion formats where ref or alt might be empty.
#' 
#' For SNV or MNV, usually no change is required.
#' For an insertion, the final REF should have a length of 1 (the anchor base),
#' and ALT should have length at least of 2
#' For a deletion, the final ALT should have a length of 1 (the anchor base),
#' and REF should have length at least of 2
#'
#' If ref or alt is empty (depending on whether it is an insertion or deletion),
#' the function will fetch the anchoring base from the fasta at position pos,
#' which is assumed to be the coordinate immediately before the indel in standard VCF conventions.
#' After ensuring one shared base (the anchor), the function trims any remaining common prefix
#' so that only the first nucleotide remains the same on both ref and alt.
#'
#' @param chr A character string specifying the chromosome name "chr1".
#' @param pos An integer specifying the position of the mutation (1-based).
#' @param ref A character string specifying the reference allele (can be empty for some formats).
#' @param alt A character string specifying the alternate allele (can be empty for some formats).
#' @param mutation_type A character string indicating the type of mutation. 
#'   Valid values include: "SNV", "MNV", "insertion", "deletion".
#' @param fasta An object or path to a FASTA file. You may need to adapt the code inside
#'   to properly retrieve the reference base (using {Biostrings::getSeq} or another method).
#'
#' @return A character vector of length 2, c(new_ref, new_alt), containing the normalized reference and alternate alleles.
#'
#' @examples 
#' Example: normalizing an insertion where ref is empty and alt is "TC"
#' Position 3 is assumed to be the base before the actual insertion.
#' Suppose the anchor base in FASTA at position 3 is "C".
#' 
#' c(new_ref, new_alt) <- normalize_ref_alt(
#'   pos = 3,
#'   ref = "",
#'   alt = "TC",
#'   mutation_type = "insertion",
#'   fasta = "path/to/genome.fasta"
#' )
#' Result might be new_ref = "C", new_alt = "CTC" (depending on trimming).
#'

normalize_ref_alt <- function(chr, pos, ref, alt, mutation_type, fasta) {
  
  # If the mutation_type is SNV or MNV, usually no normalization is required.
  # Return as-is in that case.
  if (mutation_type %in% c("SNV", "MNV")) {
    return(c(ref, alt))
  }
  
  # For insertion or deletion, we need to ensure the standard VCF convention.
  # 1) Possibly fetch anchor base from FASTA if ref or alt is empty.
  # 2) Prepend that anchor base to both ref and alt if needed.
  # 3) Trim any common prefix so that only the first nucleotide remains identical.
  
  if (mutation_type == "insertion") {
    # If REF is empty for an insertion, fetch the anchor base and prepend it.
    if (nchar(ref) == 0) {
      anchor_base <- retrieve_anchor_base(chr, pos, fasta)
      ref <- paste0(anchor_base, ref)
      alt <- paste0(anchor_base, alt)
    }
  } else if (mutation_type == "deletion") {
    # If ALT is empty for a deletion, fetch the anchor base and prepend it.
    if (nchar(alt) == 0) {
      anchor_base <- retrieve_anchor_base(chr, pos, fasta)
      ref <- paste0(anchor_base, ref)
      alt <- paste0(anchor_base, alt)
    }
  }
  
  trimmed <- trim_left_prefix(ref, alt)
  ref <- trimmed$a
  alt <- trimmed$b
  
  return(c(ref, alt))
}
