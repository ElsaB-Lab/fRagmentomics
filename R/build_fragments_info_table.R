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
build_fragments_info_table <- function(df_sam,
                                       sample_id,
                                       chr,
                                       pos,
                                       ref,
                                       alt,
                                       del_info,
                                       mutation_type,
                                       n_cores = 1) {
  
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

