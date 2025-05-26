# ------------------------------------------------------
# Helper Functions for FASTA setup and cleanup
# ------------------------------------------------------
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

# ------------------------------------------------------
# Tests
# ------------------------------------------------------
test_that("Basic insertion (AT repeat) works as expected", {
  sequences <- c(chr1 = "CATATATATATAG")
  fasta_env <- setup_test_fasta(sequences)
  on.exit(cleanup_test_fasta(fasta_env), add = TRUE)

  result <- get_repetition_seq_info(
    chr = "chr1", pos = 1, ref = "C", alt = "CAT",
    fasta = fasta_env$fa_obj, mutation_type = "insertion"
  )
  expect_equal(result, 13)
})

test_that("Basic deletion (AG repeat) works as expected", {
  sequences <- c(chr_del = "CAGAGAGAGT")
  fasta_env <- setup_test_fasta(sequences)
  on.exit(cleanup_test_fasta(fasta_env), add = TRUE)

  result <- get_repetition_seq_info(
    chr = "chr_del", pos = 1, ref = "CAG", alt = "C",
    fasta = fasta_env$fa_obj, mutation_type = "deletion"
  )
  expect_equal(result, 10)
})

test_that("Insertion with no repetition (immediate mismatch on second indel base)", {
  sequences <- c(chr_norep = "CAGTC")
  fasta_env <- setup_test_fasta(sequences)
  on.exit(cleanup_test_fasta(fasta_env), add = TRUE)

  result <- get_repetition_seq_info(
    chr = "chr_norep", pos = 1, ref = "C", alt = "CAT",
    fasta = fasta_env$fa_obj, mutation_type = "insertion"
  )
  expect_equal(result, 3)
})

test_that("Insertion with repetition until end of chromosome", {
  sequences <- c(chr_end = "CATAT")
  fasta_env <- setup_test_fasta(sequences)
  on.exit(cleanup_test_fasta(fasta_env), add = TRUE)

  # No warning is expected when reading off the end of a chromosome.
  # The function should gracefully handle the empty result from getSeq()
  # and return the position where the read was attempted.
  result <- get_repetition_seq_info(
    chr = "chr_end", pos = 1, ref = "C", alt = "CAT",
    fasta = fasta_env$fa_obj, mutation_type = "insertion"
  )

  expect_equal(result, 6)
})

test_that("SNV mutation_type returns NULL with warning", {
  sequences <- c(chr_snv = "CGT")
  fasta_env <- setup_test_fasta(sequences)
  on.exit(cleanup_test_fasta(fasta_env), add = TRUE)

  expect_warning(
    result <- get_repetition_seq_info(
      chr = "chr_snv", pos = 1, ref = "C", alt = "G",
      fasta = fasta_env$fa_obj, mutation_type = "SNV"
    ),
    regexp = "mutation_type 'SNV' is not 'insertion' or 'deletion'"
  )
  expect_null(result)
})

test_that("Empty indel_seq from insertion returns NULL with warning", {
  sequences <- c(chr_empty = "GATTACA")
  fasta_env <- setup_test_fasta(sequences)
  on.exit(cleanup_test_fasta(fasta_env), add = TRUE)

  expect_warning(
    result <- get_repetition_seq_info(
      chr = "chr_empty", pos = 1, ref = "C", alt = "A",
      fasta = fasta_env$fa_obj, mutation_type = "insertion"
    ),
    regexp = "Derived indel sequence is empty"
  )
  expect_null(result)
})

test_that("Empty indel_seq from deletion returns NULL with warning", {
  sequences <- c(chr_empty_del = "GATTACA")
  fasta_env <- setup_test_fasta(sequences)
  on.exit(cleanup_test_fasta(fasta_env), add = TRUE)

  expect_warning(
    result <- get_repetition_seq_info(
      chr = "chr_empty_del", pos = 1, ref = "A", alt = "C",
      fasta = fasta_env$fa_obj, mutation_type = "deletion"
    ),
    regexp = "Derived indel sequence is empty"
  )
  expect_null(result)
})

test_that("Chromosome not found in FASTA returns NULL with specific warning", {
  sequences <- c(chr_exists = "GATTACA")
  fasta_env <- setup_test_fasta(sequences)
  on.exit(cleanup_test_fasta(fasta_env), add = TRUE)

  expect_warning(
    result <- get_repetition_seq_info(
      chr = "chr_nonexistent",
      pos = 10,
      ref = "C",
      alt = "CAT",
      fasta = fasta_env$fa_obj,
      mutation_type = "insertion"
    ),
    regexp = "Details:.*record.*failed"
  )
  expect_null(result)
})

test_that("Function handles indel with standard VCF-like representation", {
  sequences <- c(chr_start_anchor = "CATATATAG")
  fasta_env <- setup_test_fasta(sequences)
  on.exit(cleanup_test_fasta(fasta_env), add = TRUE)

  result <- get_repetition_seq_info(
    chr = "chr_start_anchor", pos = 1, ref = "C", alt = "NAT",
    fasta = fasta_env$fa_obj, mutation_type = "insertion"
  )
  expect_equal(result, 9)
})
