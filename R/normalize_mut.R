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
  df_mut_norm <- data.frame(
    chr = character(),
    pos = integer(),
    ref = character(),
    alt = character()
  )
  rename_map <- c(chr = "CHROM", pos = "POS", ref = "REF", alt = "ALT")

  for (i in seq_len(nrow(df_mut))) {
    chr <- df_mut[i, "CHROM"]
    pos <- df_mut[i, "POS"]
    ref <- df_mut[i, "REF"]
    alt <- df_mut[i, "ALT"]

    # Extract initial mutation informations
    input_mutation_info <- paste0(chr, ":", pos, ":", ref, "-", alt)

    # Normalization user-provided representation into vcf representation
    df_mut_vcf_norm <- normalize_to_vcf_rep(
      chr = chr, pos = pos, ref = ref,
      alt = alt, fasta_fafile = fasta_fafile, one_based = one_based,
      verbose = verbose
    )

    # Sanity check to see if ref != fasta
    if (is.null(df_mut_vcf_norm)) {
      next
    }

    # Normalization vcf representation with bcftools norm
    if (apply_bcftools_norm) {
      df_mut_bcftools_norm <- apply_bcftools_norm(
        chr = df_mut_vcf_norm$chr, pos = df_mut_vcf_norm$pos,
        ref = df_mut_vcf_norm$ref, alt = df_mut_vcf_norm$alt,
        fasta = fasta, tmp_folder = tmp_folder, verbose = verbose
      )

      # Sanity check to see if bcftools worked properly
      if (is.null(df_mut_bcftools_norm)) {
        next
      }

      # rename columns
      for (i in seq_along(rename_map)) {
        old_name <- rename_map[i]
        new_name <- names(rename_map)[i]
        if (old_name %in% colnames(df_mut_bcftools_norm)) {
          colnames(df_mut_bcftools_norm)[
            colnames(df_mut_bcftools_norm) == old_name
          ] <- new_name
        }
      }

      # Add the original mutation info string as a new column
      df_mut_bcftools_norm$input_mutation_info <- input_mutation_info

      # Append to the final dataframe
      df_mut_norm <- rbind(df_mut_norm, df_mut_bcftools_norm)
    } else {
      # Add the original mutation info string as a new column
      df_mut_vcf_norm$input_mutation_info <- input_mutation_info

      # Append to the final dataframe
      df_mut_norm <- rbind(df_mut_norm, df_mut_vcf_norm)
    }
  }

  df_mut_norm
}
