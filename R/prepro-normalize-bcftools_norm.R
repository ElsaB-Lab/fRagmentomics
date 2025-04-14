# Project: ElsaBLab_fRagmentomics

#' Create a temporary VCF file
#' This function creates a temporary VCF file.
#'
#' @inheritParams apply_bcftools_norm
#'
#' @return Character. The path to the temporary VCF file.
#'
#' @noRd
create_temporary_vcf <- function(chr, pos, ref, alt) {
  # Create a temporary file with .vcf extension
  tmp_vcf <- tempfile(fileext = ".vcf")

  # Define VCF header
  vcf_header <- paste0(
    "##fileformat=VCFv4.2\n",
    "##contig=<ID=", chr, ">\n",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  )

  # Define VCF content (one variant per line)
  vcf_content <- paste(chr, pos, ".", ref, alt, ".", ".", ".", sep = "\t")

  # Write the VCF header and content to the temporary file
  writeLines(c(vcf_header, vcf_content), tmp_vcf)

  tmp_vcf
}

#' Normalize VCF file using bcftools norm
#' This function normalizes a VCF file by calling bcftools norm on a temp VCF.
#' The normalized output is then read into a dataframe, extracting the columns:
#' chr, pos, ref, and alt.
#'
#' @inheritParams process_fragment
#' @param fasta A reference genome in FASTA format.
#' @param tmp_folder Character vector for the folder temporary path.
#'
#' @return A dataframe containing normalized variant information with bcftools.
#'
#' @keywords internal
apply_bcftools_norm <- function(chr, pos, ref, alt, fasta, tmp_folder) {
  # Pass if REF and ALT are equals
  if (ref == alt) {
    normalized_variants <- data.frame(
      chr,
      pos,
      ref,
      alt,
      stringsAsFactors = FALSE
    )
  } else {
    # Create a temporary VCF file from the input variant information
    tmp_vcf <- create_temporary_vcf(chr, pos, ref, alt)

    # Create a temporary file to store the output of bcftools norm
    tmp_out_vcf <- tempfile(tmpdir = tmp_folder, fileext = ".vcf")

    # Ensure that temporary files are removed after the function executes
    on.exit(
      {
        if (file.exists(tmp_vcf)) file.remove(tmp_vcf)
        if (file.exists(tmp_out_vcf)) file.remove(tmp_out_vcf)
      },
      add = TRUE
    )

    # Build the bcftools norm command
    cmd <- sprintf(
      "bcftools norm -m +both -d exact --check REF,ALT -f %s -o %s %s",
      shQuote(fasta),
      shQuote(tmp_out_vcf),
      shQuote(tmp_vcf)
    )

    # Execute the command to normalize the VCF file
    exit_status <- tryCatch(
      {
        system(cmd, intern = FALSE, ignore.stderr = FALSE)
        0 # Success
      },
      error = function(e) {
        warning(sprintf(
          "Error running bcftools norm for %s:%d:%s:%s - %s",
          chr,
          pos,
          ref,
          alt,
          e$message
        ))
        return(1) # Fail
      }
    )

    if (exit_status != 0) {
      return(NULL)
    } else {
      # Check if output vcf is null
      if (!file.exists(tmp_out_vcf) || file.info(tmp_out_vcf)$size == 0) {
        warning(sprintf(
          "bcftools norm produced an empty VCF for %s:%d:%s:%s",
          chr,
          pos,
          ref,
          alt
        ))
        return(NULL)
      }

      # Read normalized VCF into a df
      # On capture l'éventuelle erreur si le parsing échoue
      normalized_variants <- parser_vcf(tmp_out_vcf)
    }
  }

  # Return the dataframe with normalized variant data (ou NULL si erreur)
  normalized_variants
}
