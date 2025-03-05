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
normalize_variants <- function(chr, ref, alt) {
  # Remove only the forbidden characters (-, ., _, NA) while conserving the rest of the sequence
  ref <- gsub("[-._]|NA", "", ref)
  alt <- gsub("[-._]|NA", "", alt)
  
  # Conversion into string and add chr if necessary 
  chr <- as.character(chr)
  chr[!grepl("^chr", chr)] <- paste0("chr", chr[!grepl("^chr", chr)])
  
  # Return variables
  return(list(chr = chr, ref = ref, alt = alt))
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

#' Detect chromosome naming convention in a FASTA reference
#'
#' This function checks the naming convention of chromosomes in a multi-chromosome
#' FASTA reference file. It determines whether chromosome names are prefixed with "chr"
#' (e.g., "chr1") or simply numeric (e.g., "1").
#'
#' @param fasta A character string specifying the path to a FASTA reference file.
#'
#' @return A logical value:
#'   - `TRUE` if the FASTA uses "chr"-style naming (e.g., "chr1").
#'   - `FALSE` if the FASTA uses numeric naming (e.g., "1").
#'   - `NA` if the function cannot determine the naming convention.
#' 
#' @noRd
detect_fasta_chr_style <- function(fasta) {
  
  # Fetch the names (i.e., FASTA headers) for each sequence
  idx <- scanFaIndex(fasta)
  seq_names <- as.character(seqnames(idx))

  # Return TRUE if "chr1" is found
  if ("chr1" %in% seq_names) {
    return(TRUE)
    
  # Return FALSE if "1" is found
  } else if ("1" %in% seq_names) {
    return(FALSE)
    
  # Otherwise return NA, with a warning
  } else {
    return(NA)
  }
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

  # # Convert a 0-based position to a 1-based position
  if (!one_based) {
    pos <- pos + 1
  }

  # Normalize chr column into chr1 convention 
  # Remove all the forbidden characters (-, ., _, NA) while conserving the rest of the sequence
  result <- normalize_variants(chr, ref, alt)

  # Get back the variables
  chr_norm <- result$chr
  ref_norm <- result$ref
  alt_norm <- result$alt

  # Check if chr_norm is in chr... format
  if (!grepl("^chr", chr_norm)) {
    chr_norm <- paste0("chr", chr_norm)  # Convert "1" → "chr1"
  }

  # Load fasta
  fasta_loaded <- FaFile(fasta)
  open(fasta_loaded)

  # Adapt the chr convention to search the fasta
  if (isTRUE(detect_fasta_chr_style(fasta_loaded))) {
    chr_norm_fasta <- chr_norm
  } else if (isFALSE(detect_fasta_chr_style(fasta_loaded))) {
    # FASTA uses numeric format → Ensure chr_norm also follows this format
    chr_norm_fasta <- sub("^chr", "", chr_norm) 
  } else {
    warning(paste0("Unable to determine naming convention (chr vs. numeric) in the fasta."))
    return(NULL)
  }

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
        chr      = chr_norm_fasta,
        position = pos - 1, # shift one base to the left to find the anchor
        fasta    = fasta_loaded
      )

      # Add this base before the sequence of ref and alt 
      alt_norm <- anchor_base
      ref_norm <- paste0(anchor_base, ref_norm)
      pos      <- pos -1 

    } 
  } else {
    # Insertion. In this case, we add the nucleotide at the position of interest. 
    # Indeed, position meant to match the position before the insertion in the ref
    if (nchar(ref_norm) == 0) {
        anchor_base <- retrieve_anchor_base(
        chr      = chr_norm_fasta,
        position = pos,
        fasta    = fasta_loaded
      )

      # Add this base before the sequence of ref and alt 
      alt_norm <- paste0(anchor_base, alt_norm)
      ref_norm <- anchor_base
    } 
  }

  # Sanity check: Seq ref = fasta at the position 
  check_result <- sanity_check_prepro_user_to_vcf_repr(chr_norm_fasta, pos, ref_norm, fasta_loaded)
  if (check_result) {
    message("Reference allele matches FASTA")
    # Return the final normalized values
    return(list(chr = chr_norm, pos = pos, ref = ref_norm, alt = alt_norm))
  } else {
    warning(paste0("Mismatch found between ref and fasta for (", chr, " ", pos, " ", ref, ")."))
    return(NULL)
  }
}
