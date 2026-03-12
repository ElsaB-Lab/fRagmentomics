#' Normalize a data frame of variants
#'
#' @description This function serves as the main variant normalization pipeline.
#' It iterates through a data frame of variants, applying a two-step
#' normalization process to each one, ensuring they are represented in a
#' canonical, left-aligned format suitable for downstream analysis.
#'
#' @inheritParams run_fRagmentomics
#' @param df_mut a dataframe with mutations
#' @param fasta_fafile An open connection to an object of class FaFile
#'
#' @return a dataframe with normalized mutations representations containing
#' columns CHROM, POS, REF, ALT
#'
#' @keywords internal
normalize_mut <- function(
  df_mut, fasta, fasta_fafile, one_based,
  apply_bcftools_norm, tmp_folder, verbose
) {
  rename_map <- c(chr = "CHROM", pos = "POS", ref = "REF", alt = "ALT")

  results <- lapply(seq_len(nrow(df_mut)), function(i) {
    chr <- df_mut[i, "CHROM"]
    pos <- df_mut[i, "POS"]
    ref <- df_mut[i, "REF"]
    alt <- df_mut[i, "ALT"]

    # Store original mutation information string for traceability
    input_mutation_info <- paste0(chr, ":", pos, ":", ref, "-", alt)

    # Normalize user-provided representation to VCF format
    df_mut_vcf_norm <- normalize_to_vcf_rep(
      chr = chr, pos = pos, ref = ref,
      alt = alt, fasta_fafile = fasta_fafile, one_based = one_based,
      verbose = verbose
    )

    # Sanity check to see if ref != fasta
    if (is.null(df_mut_vcf_norm)) {
      return(NULL)
    }

    # Further normalize to left-aligned VCF representation using bcftools norm
    if (apply_bcftools_norm) {
      df_mut_bcftools_norm <- apply_bcftools_norm(
        chr = df_mut_vcf_norm$chr, pos = df_mut_vcf_norm$pos,
        ref = df_mut_vcf_norm$ref, alt = df_mut_vcf_norm$alt,
        fasta = fasta, tmp_folder = tmp_folder, verbose = verbose
      )

      # Sanity check to see if bcftools worked properly
      if (is.null(df_mut_bcftools_norm)) {
        return(NULL)
      }

      # Rename columns to match the expected output format (vectorized)
      idx <- match(rename_map, colnames(df_mut_bcftools_norm))
      valid_idx <- !is.na(idx)
      colnames(df_mut_bcftools_norm)[idx[valid_idx]] <- names(rename_map)[valid_idx]

      df_mut_bcftools_norm$input_mutation_info <- input_mutation_info
      df_mut_bcftools_norm
    } else {
      df_mut_vcf_norm$input_mutation_info <- input_mutation_info
      df_mut_vcf_norm
    }
  })

  filtered <- Filter(Negate(is.null), results)
  if (length(filtered) == 0) {
    return(data.frame(
      chr = character(), pos = integer(), ref = character(),
      alt = character(), input_mutation_info = character()
    ))
  }
  do.call(rbind, filtered)
}
