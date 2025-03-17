# Project : ElsaBLab_fRagmentomics

#' Expand Multiallelic Variants
#'
#' This function expands multiallelic variants by splitting the ALT column at commas,
#' creating a new row for each allele while preserving CHROM (in string), POS, and REF values.
#'
#' @param df A data frame containing variant data with columns CHROM, POS, REF, and ALT.
#' @return A data frame where each row contains only a single alternative allele.
#' 
#' @noRd
expand_multiallelics <- function(df) {
  expanded_rows <- list()

  for (i in seq_len(nrow(df))) {
    # Before doing the strsplit, we have to be sure ALT doesn't finish by ","
    if (grepl(",$", df$ALT[i])) {
        df$ALT[i] <- paste0(df$ALT[i], "-")
    }
    
    # Devide all the alternatives alleles
    alts <- unlist(strsplit(df$ALT[i], ","))  

    for (alt in alts) {
      expanded_rows <- append(expanded_rows, list(
        data.frame(CHROM = df$CHROM[i], POS = df$POS[i], REF = df$REF[i], ALT = alt, stringsAsFactors = FALSE)
      ))
    }
  }

  # Recreate the df with one line for each alt 
  expanded_df <- do.call(rbind, expanded_rows)
  return(expanded_df)
}

#' Read VCF File
#'
#' This function reads a VCF file, applies sanity checks, and expands multiallelic variants.
#'
#' @param vcf_file Path to the VCF file.
#' @return A processed data frame with normalized variants.
#' 
#' @noRd
read_vcf_input <- function(vcf_file) {
    vcf_subset <- parser_vcf(vcf_file)
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
    tsv_subset <- parser_tsv(tsv_file)
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
read_str_input <- function(mut) {
    mut_df <- parser_mut_str(mut)
    mut_df_without_multiallelic <- expand_multiallelics(mut_df)
    return(mut_df_without_multiallelic)
}


#' Get the informations about mutations
#'
#' @inheritParams fRagmentomics
#' @param mut could be a .vcf, .tsv or chr:pos:ref:alt
#' 
#' @return a df without multiallelic ALT and without header with chr pos ref alt columns
#'
#' @noRd
read_mut <- function(mut) {
  if (grepl("\\.tsv(.gz)?$", mut)) {

    mut_df <- read_tsv_input(mut)

  } else if (grepl("\\.vcf(.gz)?$", mut)) {
    
    mut_df <- read_vcf_input(mut)

  } else if (grepl("^[^:]+:[^:]+:[^:]+:[^:]+$", mut)) {
    
    mut_df <- read_str_input(mut)

  } else {
    stop(paste0("Error: The parameter 'mut' (", mut, ") is not in the expected format (.tsv, .vcf, chr:pos:ref:alt)."))
  }
  
  return(mut_df)
}


