# Project : ElsaBLab_fRagmentomics

#' @title Setup parallel computations
#' @description Initializes and registers a parallel computing cluster for processing.
#'
#' @param n_cores Number of cores to use for parallel processing.
#' @return A parallel cluster object.
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' 
#' @noRd
setup_parallel_computations <- function(n_cores){
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  return(cl)
}

#' @title Build fragment information table
#' @description Processes fragment-level sequencing data in parallel and extracts relevant statistics.
#'
#' @param df_sam A dataframe containing sequencing reads.
#' @param sample_id Sample identifier.
#' @param chr Chromosome of interest.
#' @param pos Genomic position of interest.
#' @param ref Reference base.
#' @param alt Alternative base (for SNV) or sequence (for insertion).
#' @param del_info Deletion-specific information (if applicable).
#' @param mutation_type Type of mutation: "mutation", "deletion", or "insertion".
#' @param n_cores Number of cores for parallel computation.
#'
#' @return A dataframe containing extracted fragment-level information.
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel stopCluster
#' 
#' @export
build_fragments_info_table <- function(
    bam,
    mut,
    fasta,
    sample = NA,
    neg_offset_mate_search = -1000,
    pos_offset_mate_search = 1000,
    flag_keep = 0x03,
    flag_remove = 0x900, 
    report_5p_bases = 0,
    report_3p_bases = 0,
    report_tlen = FALSE,
    report_softclip = FALSE, 
    n_cores = 1
  ) {
  
  # Read a vcf ou .tsv file 
  # Return a df with all the mutation we want to study 
  # Look at the format of mut to know if it is a VCF, mut_file, chr:pos:ref:alt
  mut_info <- read_mut(mut)

  # Loop on each row of the mut_info
  for (i in 1:nrow(mut_info)) {
    chr <- mut_info[i, 1]
    pos <- mut_info[i, 2]
    ref <- mut_info[i, 3]
    alt <- mut_info[i, 4]

  # Read and extract bam around the mutation position --> Boucle  
  # Return a troncated sam 
  sam <- preprocess_bam(
    bam,
    chr,
    pos,
    neg_offset_mate_search,
    pos_offset_mate_search,
    flag_keep,
    flag_remove)

  # Find mutation_status "insertion", "deletion", "SNV", "MNP"
  mutation_status <- define_mutation_status(ref, alt)

  # Normalize the REF and ALT format for the indel. We have the mutation status 
  c(ref, alt) <- normalize_ref_alt(ref, alt)

  # Process fragmentomics on truncated bam and the mutation 
  process_fragmentomics <- function(
    sam,
    sample_id,
    chr,
    pos,
    ref,
    alt,
    out,
    n_cores = 1
  ) {




  # Extract unique fragment names
  fragments_names <- unique(df_sam[, 1, drop = TRUE])
  n_fragments <- length(fragments_names)
  
  # Initialize parallel cluster
  cl <- setup_parallel_computations(n_cores)
  
  # Parallel execution
  df_fragments_info <- foreach::foreach(
    i = seq_len(n_fragments),
    .combine = rbind,
    .inorder = FALSE,
    .packages = c("stringr")
  ) %dopar% {
    # Call the fragment processing function
    process_fragment(
      df_sam           = df_sam,
      fragment_name    = fragments_names[i],
      sample_id        = sample_id,
      chr              = chr,
      pos              = pos,
      ref              = ref,
      alt              = alt,
      del_info         = del_info,
      mutation_type    = mutation_type
    )
  }
  
  # Stop cluster
  parallel::stopCluster(cl)
  
  return(df_fragments_info)
  }
  } 
}

