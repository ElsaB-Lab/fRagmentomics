# Project : ElsaBLab_fRagmentomics

#' Read VCF File
#'
#' This function reads a VCF file, applies sanity checks, and expands multiallelic variants.
#'
#' @param vcf_file Path to the VCF file.
#' @return A processed data frame with normalized variants.
#' 
#' @noRd
read_vcf_input <- function(vcf_file) {
    vcf_subset <- sanity_check_vcf(vcf_file)
    # Create a row for each alt in multiallelics cases
    vcf_subset_without_multiallelic <- expand_multiallelics(vcf_subset)
    return(vcf_subset_without_multiallelic)
}

#' Read TSV File
#'
#' This function reads a TSV file, applies sanity checks, and expands multiallelic variants.
#'
#' @param tsv_file Path to the TSV file.
#' @return A processed data frame with normalized variants.
#' 
#' @noRd
read_tsv_input <- function(tsv_file) {
    tsv_subset <- sanity_check_tsv(tsv_file)
    tsv_subset_without_multiallelic <- expand_multiallelics(tsv_subset)
    return(tsv_subset_without_multiallelic)
}

#' Read mut informations
#'
#' This function reads a chr:pos:ref:alt format
#'
#' @param mut Infos chr:pos:ref:alt
#' @return A processed data frame with normalized variants.
#' 
#' @noRd
read_mut_input <- function(mut) {
    # If mut is finished by ":", the programm will add - at the end to use strsplit
    if (grepl(":$", mut)) {
      mut <- paste0(mut, "-")
    }

    # Devide the 4 parts of the mutation expected
    parts <- unlist(strsplit(mut, ":"))

    # Check if mutation contains exactly these parts (CHR, POS, REF, ALT)
    if (length(parts) != 4) {
      stop(paste0("Error: The parameter 'mut' (", mut, ") is not in the expected format (.tsv, .vcf, chr:pos:ref:alt)."))
    }

    chr <- parts[1]
    pos <- pos <- as.integer(parts[2])
    ref <- ref <- parts[3]
    alt <- alt <- parts[4]

    # Return the info into a list 
    mut_df <- data.frame(CHROM = chr, POS = pos, REF = ref, ALT = alt, stringsAsFactors = FALSE)

    # We take into account multiallelic in ALT column 
    mut_df_without_multiallelic <- expand_multiallelics(mut_df)
    return(mut_df_without_multiallelic)
}


#' Get the informations about mutations
#'
#' @inheritParams fRagmentomics
#' @param mut could be a .vcf, .tsv or chr:pos:ref:alt
#' @return a df without multiallelic ALT
#'
#' @noRd
read_mut <- function(mut) {
  # Define the REGEX that will be used to capture chr:pos:ref:alt
  chromosome_pattern <- "([0-9XY]+|chr[0-9XY]*)"             # Mandatory: chr1 or 1
  position_pattern <- ":[0-9]+"                              # Mandatory: A number
  ref_pattern <- "([ACGT._-]+|NA)"                           # REF can be ACGT, ".", "_" o "-"
  alt_pattern <- "([ACGT._-]+(,[ACGT._-]*)*|NA)"             # ALT can be ACGT, ".", "_" o "-"

  # Different possible cases
  ref_alt_pattern <- paste0(ref_pattern, ":", alt_pattern)   # Normal case (REF:ALT)
  ref_missing_pattern <- paste0(ref_pattern, ":")            # REF but no ALT (REF:)
  alt_missing_pattern <- paste0(":", alt_pattern)            # ALT but no REF (:ALT)

  # Combine all patterns into one full regex
  full_pattern <- paste0("^", chromosome_pattern, position_pattern, 
                        ":(?:", ref_alt_pattern, "|", ref_missing_pattern, "|", alt_missing_pattern, ")$")

  # Return a datafram without header, with all the mutations from the vcf
  # Columns : chr pos ref alt

  if (grepl("\\.tsv$", mut)) {

    mut_df <- read_tsv_input(mut)

  } else if (grepl("\\.vcf$", mut)) {
    
    mut_df <- read_vcf_input(mut)

  } else if (grepl(full_pattern, mut)) {
    
    mut_df <- read_mut_input(mut)

  } else {
    stop(paste0("Error: The parameter 'mut' (", mut, ") is not in the expected format (.tsv, .vcf, chr:pos:ref:alt)."))
  }
  
  return(mut_df)
}


