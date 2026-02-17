#' Check that bcftools system dependency is available.
#'
#' @return The system path to the bcftools executable
#' @keywords internal
check_bcftools_is_installed <- function() {
  bcftools_path <- Sys.which("bcftools")
  if (bcftools_path == "") {
    stoptext <- paste0(
      "bcftools is not in your PATH.\nPlease install it and ensure",
      "it's available in your system's PATH.\n",
      "See the package vignette for installation instructions."
    )
    stop(stoptext)
  }
  return(bcftools_path)
}


#' Normalize a single variant using bcftools norm
#'
#' @description This function normalizes a single variant by leveraging the
#' external 'bcftools norm' command. It writes the variant to a temporary VCF
#' file, executes 'bcftools norm' for left-alignment and parsimonious
#' representation, and then reads the normalized result back into a data frame.
#'
#' @inheritParams normalize_mut
#' @param chr A string representing the chromosome.
#' @param pos An integer representing the position.
#' @param ref A string representing the reference allele.
#' @param alt A string representing the alternative allele.
#'
#' @return A dataframe containing normalized variant information with bcftools.
#'
#' @keywords internal
apply_bcftools_norm <- function(
  chr, pos, ref, alt, fasta, tmp_folder,
  verbose
) {
  # Pass if REF and ALT are equals
  if (ref == alt) {
    normalized_variants <- data.frame(chr, pos, ref, alt,
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
    command <- "bcftools"
    args <- c(
      "norm", "-m", "+both", "-d", "exact", "--check", "REF,ALT",
      "-f", fasta, "-o", tmp_out_vcf, tmp_vcf
    )

    # Print a message
    if (verbose) {
      messtext <- sprintf(
        paste0(
          "Executing bcftools norm normalisation for variant",
          " for %s:%d:%s:%s"
        ), chr, pos, ref, alt
      )
      message(messtext)
    }

    # Execute the command to normalize the VCF file
    exit_status <- tryCatch(
      {
        # Execute the command
        system2(command, args = args, stdout = TRUE, stderr = TRUE)
        0 # Success
      },
      error = function(e) {
        warning(sprintf(
          "Bcftools norm failed for variant %s:%d:%s:%s - %s",
          chr, pos, ref, alt, e$message
        ))
        return(1) # Fail
      }
    )

    if (exit_status != 0) {
      return(NULL)
    } else {
      # Check if output vcf is null
      if (!file.exists(tmp_out_vcf) || file.info(tmp_out_vcf)$size == 0) {
        warntext <- sprintf(
          paste0(
            "bcftools norm produced an empty VCF",
            "for %s:%d:%s:%s"
          ), chr, pos, ref, alt
        )
        warning(warntext)
        return(NULL)
      }

      # Read normalized VCF into a df
      normalized_variants <- parser_vcf(tmp_out_vcf)
    }
  }

  # Return the dataframe with normalized variant data (or NULL if erreur)
  normalized_variants
}


#' Create a temporary VCF file
#'
#' @description This function creates a temporary VCF file.
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
    "##fileformat=VCFv4.2\n", "##contig=<ID=", chr, ">\n",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  )

  # Define VCF content (one variant per line)
  vcf_content <- paste(chr, pos, ".", ref, alt, ".", ".", ".", sep = "\t")

  # Write the VCF header and content to the temporary file
  writeLines(c(vcf_header, vcf_content), tmp_vcf)

  tmp_vcf
}
