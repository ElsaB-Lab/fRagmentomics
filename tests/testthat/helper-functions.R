setup_test_fasta <- function(named_sequences, dir = tempfile("testfasta_dir")) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  fasta_path <- file.path(dir, "test_genome.fasta")

  dna_strings <- Biostrings::DNAStringSet(named_sequences)
  Biostrings::writeXStringSet(dna_strings, fasta_path)
  Rsamtools::indexFa(fasta_path)

  fa_file_obj <- Rsamtools::FaFile(fasta_path)
  open(fa_file_obj)
  return(list(fa_obj = fa_file_obj, path = fasta_path, dir = dir))
}

cleanup_test_fasta <- function(fasta_setup_list) {
  if (!is.null(fasta_setup_list$fa_obj) && inherits(fasta_setup_list$fa_obj, "FaFile")) {
    try(close(fasta_setup_list$fa_obj), silent = TRUE)
  }
  unlink(fasta_setup_list$path, force = TRUE)
  unlink(paste0(fasta_setup_list$path, ".fai"), force = TRUE)
  if (dir.exists(fasta_setup_list$dir)) {
    unlink(fasta_setup_list$dir, recursive = TRUE, force = TRUE)
  }
}

setup_test_bam <- function(df_reads, ref_seq_lengths, dir = tempfile("testbam_dir")) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }

  sam_path <- file.path(dir, "test.sam")
  bam_root <- file.path(dir, "test") # Rsamtools::asBam will add the .bam extension

  # 1. Create SAM Header
  header_lines <- c("@HD\tVN:1.0\tSO:coordinate")
  for (seq_name in names(ref_seq_lengths)) {
    header_lines <- c(header_lines, sprintf("@SQ\tSN:%s\tLN:%d", seq_name, ref_seq_lengths[[seq_name]]))
  }

  # Write Header
  writeLines(header_lines, sam_path)

  # 2. Prepare SAM Body
  if (!"MAPQ" %in% names(df_reads)) df_reads$MAPQ <- 60
  if (!"RNEXT" %in% names(df_reads)) df_reads$RNEXT <- "="
  if (!"PNEXT" %in% names(df_reads)) df_reads$PNEXT <- 0
  if (!"TLEN" %in% names(df_reads)) df_reads$TLEN <- 0

  # Ensure columns are ordered correctly for SAM format standard
  cols_ordered <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL")

  # Check if all required columns are present
  df_to_write <- df_reads[, cols_ordered]

  # 3. Write Body
  write.table(
    df_to_write,
    file = sam_path,
    append = TRUE,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )

  # 4. Convert to BAM using Rsamtools
  # This automatically creates the .bam.bai index file
  bam_path <- Rsamtools::asBam(
    file = sam_path,
    destination = bam_root,
    overwrite = TRUE,
    indexDestination = TRUE
  )

  return(list(path = bam_path, dir = dir))
}

cleanup_test_bam <- function(bam_setup_list) {
  if (dir.exists(bam_setup_list$dir)) {
    unlink(bam_setup_list$dir, recursive = TRUE, force = TRUE)
  }
}

calculate_read_length <- function(cigar) {
  ops <- stringr::str_extract_all(cigar, "[[:digit:]]+(?=[MIS=X])")[[1]]
  sum(as.numeric(ops))
}

skip_if_no_bcftools <- function() {
  if (Sys.which("bcftools") == "") {
    testthat::skip("bcftools not available for testing")
  }
}
