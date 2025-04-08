# Project : ElsaBLab_fRagmentomics

#' Normalize Variant Data
#' This function normalizes the REF and ALT columns by
#' replacing invalid values with an empty string.
#'
#' @inheritParams normalize_user_rep_to_vcf_rep
#'
#' @return A normalized data frame.
#'
#' @noRd
normalize_na_representation <- function(ref, alt) {
  # Remove only the forbidden characters (-, ., _, NA)
  ref <- gsub("[-._]|NA", "", ref)
  alt <- gsub("[-._]|NA", "", alt)

  # Return variables
  list(ref = ref, alt = alt)
}

#' Retrieve a single nucleotide from a multi-chromosome FASTA
#'
#' @inheritParams normalize_user_rep_to_vcf_rep
#'
#' @return A character string representing the single nucleotide at
#' the specified chromosome and position.
#'
#' @importFrom Biostrings getSeq
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#'
#' @noRd
retrieve_anchor_base <- function(chr, position, fasta) {
  # Fetch a single nucleotide from the FASTA
  # using the chromosome, and start = end = position
  anchor_base <- Biostrings::getSeq(
    x = fasta,
    param = GenomicRanges::GRanges(
      seqnames = chr,
      ranges = IRanges::IRanges(start = position, end = position)
    )
  )

  # Transform fasta_seq into string
  anchor_base_char <- unname(as.character(anchor_base))
  return(as.character(anchor_base_char))
}

#' Harmonize Chromosome Notation with FASTA Reference
#' This function ensures that a given chr notation (e.g., `"1"` or `"chr1"`)
#' matches the convention used in the provided FASTA.
#'
#' @inheritParams normalize_user_rep_to_vcf_rep
#'
#' @importFrom Rsamtools scanFaIndex
#' @importFrom GenomeInfoDb seqnames
#'
#' @return String with the chr notation harmonized to match the FASTA reference.
#'
#' @noRd
harmonize_chr_to_fasta <- function(chr, fasta) {
  # Extract chromosome names from the FASTA index
  fasta_index <- scanFaIndex(fasta) # Get the indexed FASTA sequence info
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

#' Normalize user-provided representation into VCF representation
#' This function normalizes the REF and ALT columns into the VCF format.
#' It also checks if all nucleotides in the REF column match those from
#' a given FASTA reference.
#'
#' @inheritParams process_fragment
#' @param fasta A reference genome in FASTA format.
#' @param one_based Boolean. TRUE if fasta is in one based. False if in 0 based.
#'
#' @return A list with chr, pos, ref and alt
#'
#' @keywords internal
normalize_user_rep_to_vcf_rep <- function(
    chr,
    pos,
    ref,
    alt,
    fasta,
    one_based) {
  # Remove all the forbidden characters (-, ., _, NA)
  result <- normalize_na_representation(ref, alt)

  # # Convert a 0-based position to a 1-based position
  if (!one_based) {
    pos_norm <- pos + 1
  } else {
    pos_norm <- pos
  }

  # Conversion into string and add chr if necessary
  chr <- as.character(chr)

  # Get back the variables
  ref_norm <- result$ref
  alt_norm <- result$alt

  # Normalize the chr convention from user with the one from the fasta
  chr_norm <- harmonize_chr_to_fasta(chr, fasta)

  # Adapt ref and alt for the vcf convention in function of the mutation type
  if (nchar(ref_norm) == nchar(alt_norm)) {
    # SNV or MNP: No change to position or sequence
    # We assume the position is correct, so just return the normalized strings.
  } else if (nchar(ref_norm) > nchar(alt_norm)) {
    # Deletion: In this case, we add the nucleotide before the pos of interest.
    # Indeed, position meant to match the first position of the del in the ref
    # We also have to return the position before
    if (nchar(alt_norm) == 0) {
      anchor_base <- retrieve_anchor_base(
        chr      = chr_norm,
        position = pos_norm - 1, # shift one base to the left to find the anchor
        fasta    = fasta
      )

      # Add this base before the sequence of ref and alt
      alt_norm <- anchor_base
      ref_norm <- paste0(anchor_base, ref_norm)
      pos_norm <- pos_norm - 1
    }
  } else {
    # Insertion. In this case, we add the nucleotide at the pos of interest.
    # Indeed, position meant to match the position before the ins in the ref
    if (nchar(ref_norm) == 0) {
      anchor_base <- retrieve_anchor_base(
        chr      = chr_norm,
        position = pos_norm,
        fasta    = fasta
      )

      # Add this base before the sequence of ref and alt
      alt_norm <- paste0(anchor_base, alt_norm)
      ref_norm <- anchor_base
    }
  }

  # Sanity check: Seq ref = fasta at the position
  check_result <- sanity_check_ref_in_fasta(chr_norm, pos_norm, ref_norm, fasta)
  if (check_result) {
    message("Reference allele matches FASTA")
    # Return the final normalized values
    return(list(chr = chr_norm, pos = pos_norm, ref = ref_norm, alt = alt_norm))
  } else {
    warning(paste0("Mismatch found between ref and fasta
    for (", chr, " ", pos, " ", ref, ")."))
    return(NULL)
  }
}
