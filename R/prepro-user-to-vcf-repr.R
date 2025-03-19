# Project : ElsaBLab_fRagmentomics

#' Normalize Variant Data
#'
#' This function normalizes the REF and ALT columns by replacing invalid values with an empty string.
#' It also ensures that the CHROM column is treated as a string and in the format chr1.
#'
#' @param chr String with chromosome of interest.
#' @param ref String with reference base.
#' @param alt String with alternative base (for SNV) or sequence (for insertion).
#' @return A normalized data frame.
#' 
#' @noRd
normalize_na_representation <- function(ref, alt) {
  # Remove only the forbidden characters (-, ., _, NA) while conserving the rest of the sequence
  ref <- gsub("[-._]|NA", "", ref)
  alt <- gsub("[-._]|NA", "", alt)
  
  # Return variables
  return(list(ref = ref, alt = alt))
}

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
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' 
#' @noRd
retrieve_anchor_base <- function(chr, position, fasta) {
  # Fetch a single nucleotide from the FASTA
  # using the chromosome, and start = end = position
  anchor_base <- Biostrings::getSeq(
    x     = fasta,
    param = GRanges(seqnames = chr,
                    ranges = IRanges(start = position, end = position))
  )

  # Transform fasta_seq into string
  anchor_base_char <- unname(as.character(anchor_base))
  return(as.character(anchor_base_char))
}

#' Harmonize Chromosome Notation with FASTA Reference
#'
#' This function ensures that a given chromosome notation (e.g., `"1"` or `"chr1"`) 
#' matches the convention used in the provided FASTA reference file.
#'
#' @param chr A character string representing a chromosome (e.g., `"1"`, `"chr1"`).
#' @param fasta A character string specifying the path to the FASTA reference file.
#'
#' @importFrom Rsamtools scanFaIndex
#' 
#' @return A character string with the chromosome notation harmonized to match the FASTA reference.
#' 
#' @noRd 
harmonize_chr_to_fasta <- function(chr, fasta) {
  # Extract chromosome names from the FASTA index
  fasta_index <- scanFaIndex(fasta)  # Get the indexed FASTA sequence info
  fasta_chromosomes <- seqnames(fasta_index)  # Extract chromosome names

  # Determine whether FASTA uses "chr1" or "1"
  fasta_format <- ifelse(any(grepl("^chr", fasta_chromosomes)), "chr", "no_chr")
  
  # Harmonize chromosome notation
  if (fasta_format == "chr" && !grepl("^chr", chr)) {
    chr <- paste0("chr", chr)  # Add "chr" if missing
  } else if (fasta_format == "no_chr" && grepl("^chr", chr)) {
    chr <- sub("^chr", "", chr)  # Remove "chr" if present
  }
  
  return(chr)
}

#' Normalize user-provided representation into VCF representation 
#'
#' This function normalizes the REF and ALT columns into the VCF format.
#' It also checks if all nucleotides in the REF column match those from a given FASTA reference.
#'
#' @param chr A string representing the chromosome (e.g., "1" or "chr1").
#' @param pos An integer representing the genomic position of the variant.
#' @param ref A string representing the reference allele at the given position.
#' @param alt A string representing the alternative allele (single, not multiple).
#' @param fasta A reference genome in FASTA format, used to validate the REF nucleotide.
#' @param one_based A logical indicating whether the input position is 1-based (TRUE) or 0-based (FALSE).
#'
#' @return A list with chr, pos, ref and alt 
#'
#' @noRd
normalize_user_rep_to_vcf_rep <- function(chr, pos, ref, alt, fasta, one_based) {
  # Normalize chr column into chr1 convention 
  # Remove all the forbidden characters (-, ., _, NA) while conserving the rest of the sequence
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
    # Deletion: In this case, we add the nucleotide before the position of interest. 
    # Indeed, position meant to match the first position of the deletion in the ref
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
      pos_norm <- pos_norm -1 

    } 
  } else {
    # Insertion. In this case, we add the nucleotide at the position of interest. 
    # Indeed, position meant to match the position before the insertion in the ref
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
    warning(paste0("Mismatch found between ref and fasta for (", chr, " ", pos, " ", ref, ")."))
    return(NULL)
  }
}
