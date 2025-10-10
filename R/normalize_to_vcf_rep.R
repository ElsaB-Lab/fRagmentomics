#' Normalize a variant to VCF representation
#'
#' @description
#' This function converts a variant from a user-provided format into the
#' standard VCF representation, which is crucial for handling indels. It also
#' harmonizes chromosome names and validates the reference allele against a
#' FASTA file.
#'
#' @details
#' The normalization process involves several key steps:
#' 1.  **Coordinate Adjustment**: Converts the input 'pos' from 0-based to 1-based if 'one_based' is 'FALSE'.
#' 2.  **Chromosome Naming**: Standardizes the chromosome name (e.g., '1' vs 'chr1') to match the provided FASTA
#'     reference using 'harmonize_chr_to_fasta'.
#' 3.  **Indel Padding**: For insertions and deletions, it prepends an "anchor" base from the reference genome to both
#'     'ref' and 'alt' alleles to create a valid VCF-style representation. The 'pos' is adjusted accordingly for
#'     deletions. SNVs and MNVs are not modified.
#' 4.  **Validation**: After normalization, it confirms that the final 'ref' allele matches the sequence in the FASTA file at the new coordinates.
#'
#' @inheritParams normalize_mut
#' @param chr A string representing the chromosome.
#' @param pos An integer representing the position.
#' @param ref A string representing the reference allele.
#' @param alt A string representing the alternative allele.
#'
#' @return A list with chr, pos, ref and alt
#'
#' @keywords internal
normalize_to_vcf_rep <- function(
    chr,
    pos,
    ref,
    alt,
    fasta_fafile,
    one_based,
    verbose) {
  # Remove all the forbidden characters (-, ., _, NA)
  result <- normalize_na_representation(ref, alt)
  ref_norm <- result$ref
  alt_norm <- result$alt

  # Convert a 0-based position to a 1-based position
  if (!one_based) {
    pos_norm <- pos + 1
  } else {
    pos_norm <- pos
  }

  # Normalize the chr convention from user with the one from the fasta
  chr <- as.character(chr)
  chr_norm <- harmonize_chr_to_fasta(chr, fasta_fafile)

  # Adapt ref and alt for the vcf convention in function of the mutation type
  # SNV or MNV: No change to position or sequence
  # We assume the position is correct, so just return the normalized strings.
  if (nchar(ref_norm) != nchar(alt_norm)) {
    if (nchar(ref_norm) > nchar(alt_norm)) {
      # Deletion: In this case, we add the nucleotide before the pos of interest.
      # Indeed, position meant to match the first position of the del in the ref
      # We also have to return the position before
      if (nchar(alt_norm) == 0) {
        anchor_base <- get_seq_from_fasta(
          chr          = chr_norm,
          start        = pos_norm - 1, # shift one base to the left to find the anchor
          end          = pos_norm - 1, # shift one base to the left to find the anchor
          fasta_fafile = fasta_fafile
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
        anchor_base <- get_seq_from_fasta(
          chr          = chr_norm,
          start        = pos_norm,
          end          = pos_norm,
          fasta_fafile = fasta_fafile
        )

        # Add this base before the sequence of ref and alt
        alt_norm <- paste0(anchor_base, alt_norm)
        ref_norm <- anchor_base
      }
    }
  }

  # Sanity check
  ref_matches_fasta <- check_if_ref_matches_fasta(
    chr          = chr_norm,
    pos          = pos_norm,
    ref          = ref_norm,
    fasta_fafile = fasta_fafile
  )

  if (ref_matches_fasta) {
    if (verbose) {
      message("Reference allele matches FASTA")
    }
    # Return the final normalized values
    return(data.frame(chr = chr_norm, pos = pos_norm, ref = ref_norm, alt = alt_norm))
  } else {
    warning(sprintf("Mismatch found between ref and fasta for (%s %d %s).", chr, pos, ref))
    return(NULL)
  }
}


#' Clean allele strings by removing ambiguous characters
#'
#' @description A helper function that removes common representations for missing or null alleles (e.g., '-', '.', '_', 'NA')
#' from reference and alternate allele strings, converting them to empty strings ('""').
#'
#' @inheritParams normalize_to_vcf_rep
#'
#' @return A normalized dataframe.
#'
#' @noRd
normalize_na_representation <- function(ref, alt) {
  # Remove only the forbidden characters (-, ., _, NA)
  ref <- gsub("[-._]|NA", "", ref)
  alt <- gsub("[-._]|NA", "", alt)

  # Return variables
  list(ref = ref, alt = alt)
}
