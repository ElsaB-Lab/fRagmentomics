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

calculate_read_length <- function(cigar) {
    ops <- stringr::str_extract_all(cigar, "[[:digit:]]+(?=[MIS=X])")[[1]]
    sum(as.numeric(ops))
}

skip_if_no_bcftools <- function() {
  if (Sys.which("bcftools") == "") {
    testthat::skip("bcftools not available for testing")
  }
}
