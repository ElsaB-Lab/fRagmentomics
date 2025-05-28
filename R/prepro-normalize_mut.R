#' Check input files for fRagmentomics function
#'
#' @inheritParams analyze_fragments
#' @param df_mut a dataframe with mutations
#' @param fasta_fafile An open connection to an object of class FaFile
#'
#' @return a dataframe with normalized mutations representations containing columns CHROM, POS, REF, ALT
#'
#' @keywords internal
normalize_mut <- function(df_mut, fasta, fasta_fafile, one_based){
  df_mut_norm <- data.frame()

  for (i in seq_len(nrow(df_mut))) {
    chr <- df_mut[i, 1]
    pos <- df_mut[i, 2]
    ref <- df_mut[i, 3]
    alt <- df_mut[i, 4]

    # Normalization user-provided representation into vcf representation
    mut_vcf_norm <- normalize_to_vcf_rep(
      chr = chr,
      pos = pos,
      ref = ref,
      alt = alt,
      fasta_fafile = fasta_fafile,
      one_based = one_based
    )

    # Sanity check to see if ref != fasta
    if (is.null(mut_vcf_norm)) {
      next
    }

    chr_norm <- mut_vcf_norm$chr
    pos_norm <- mut_vcf_norm$pos
    ref_norm <- mut_vcf_norm$ref
    alt_norm <- mut_vcf_norm$alt

    # Normalization vcf representation with bcftools norm
    mut_bcftools_norm <- apply_bcftools_norm(
      chr        = chr_norm,
      pos        = pos_norm,
      ref        = ref_norm,
      alt        = alt_norm,
      fasta      = fasta,
      tmp_folder = tmp_folder
    )

    # Sanity check to see if bcftools worked properly
    if (is.null(mut_bcftools_norm)) {
      next
    }

    # Append to the final dataframe
    df_mut_norm <- rbind(df_mut_norm, mut_bcftools_norm)
  }

  df_mut_norm
}
