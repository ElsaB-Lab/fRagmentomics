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
  Rsamtools::indexFa(fasta_path) # Create index

  fa_file_obj <- Rsamtools::FaFile(fasta_path)
  open(fa_file_obj) # Use generic open(), relies on Rsamtools library being loaded
  return(list(fa_obj = fa_file_obj, path = fasta_path, dir = dir))
}

cleanup_test_fasta <- function(fasta_setup_list) {
  # Check if fa_obj exists and is an FaFile, and then if it's open
  if (!is.null(fasta_setup_list$fa_obj) && inherits(fasta_setup_list$fa_obj, "FaFile")) {
    if (isOpen(fasta_setup_list$fa_obj)) {
      close(fasta_setup_list$fa_obj)
    }
  }
  # Delete created files
  unlink(fasta_setup_list$path, force = TRUE)
  unlink(paste0(fasta_setup_list$path, ".fai"), force = TRUE)
  # Clean up directory
  if (dir.exists(fasta_setup_list$dir)) {
    files_in_dir <- list.files(fasta_setup_list$dir, include.dirs = TRUE, no.. = TRUE)
    if (length(files_in_dir) == 0) {
      unlink(fasta_setup_list$dir, recursive = TRUE, force = TRUE)
    }
  }
}

# ------------------------------------------------------
# Tests
# ------------------------------------------------------
test_that("Basic insertion (AT repeat) works as expected", {
  sequences <- c(chr1 = "CATATATATATAG") # C (pos 1) ATATATATATA (pos 2-12) G (pos 13)
  fasta_env <- setup_test_fasta(sequences)
  on.exit(cleanup_test_fasta(fasta_env), add = TRUE)

  result <- get_repetition_seq_info(
    chr = "chr1", pos = 1, ref = "C", alt = "CAT", # Indel sequence "AT"
    fasta = fasta_env$fa_obj, mutation_type = "insertion"
  )
  expect_equal(result, 13) # Mismatch with G at position 13
})

test_that("Basic deletion (AG repeat) works as expected", {
  sequences <- c(chr_del = "CAGAGAGAGX") # C(1) AGAGAGAG(2-9) X(10)
  fasta_env <- setup_test_fasta(sequences)
  on.exit(cleanup_test_fasta(fasta_env), add = TRUE)

  result <- get_repetition_seq_info(
    chr = "chr_del", pos = 1, ref = "CAG", alt = "C", # Indel sequence "AG"
    fasta = fasta_env$fa_obj, mutation_type = "deletion"
  )
  expect_equal(result, 10)
})

test_that("Insertion with no repetition (immediate mismatch on second indel base)", {
  sequences <- c(chr_norep = "CAGTC") # C(1) A(2) G(3) T(4) C(5)
  fasta_env <- setup_test_fasta(sequences)
  on.exit(cleanup_test_fasta(fasta_env), add = TRUE)

  # Indel "AT". Genome A matches indel A. Genome G mismatches indel T.
  result <- get_repetition_seq_info(
    chr = "chr_norep", pos = 1, ref = "C", alt = "CAT", # Indel "AT"
    fasta = fasta_env$fa_obj, mutation_type = "insertion"
  )
  expect_equal(result, 3) # Mismatch at genome position 3 (G vs T)
})

test_that("Insertion with repetition until end of chromosome", {
  sequences <- c(chr_end = "CATAT") # C(1) A(2) T(3) A(4) T(5)
  fasta_env <- setup_test_fasta(sequences)
  on.exit(cleanup_test_fasta(fasta_env), add = TRUE)

  expect_warning(
    result <- get_repetition_seq_info(
      chr = "chr_end", pos = 1, ref = "C", alt = "CAT",
      fasta = fasta_env$fa_obj, mutation_type = "insertion"
    ),
    # Match the warning from your function for a SeqReadError
    regexp = "Error trying to retrieve sequence from FASTA at chr_end:6[.] Details:.*out-of-bounds"
  )
  expect_equal(result, 6) # Failure occurs when trying to access pos 6
})

test_that("SNV mutation_type returns NULL with warning", {
  sequences <- c(chr_snv = "CGT")
  fasta_env <- setup_test_fasta(sequences)
  on.exit(cleanup_test_fasta(fasta_env), add = TRUE)

  expect_warning(
    result <- get_repetition_seq_info(
      chr = "chr_snv", pos = 1, ref = "C", alt = "G",
      fasta = fasta_env$fa_obj, mutation_type = "SNV" # Not "insertion" or "deletion"
    ),
    regexp = "do not represent an indel as they have the same length"
  )
  expect_null(result)
})

test_that("Empty indel_seq from insertion (alt allele too short) returns NULL with warning", {
  sequences <- c(chr_empty = "GATTACA")
  fasta_env <- setup_test_fasta(sequences)
  on.exit(cleanup_test_fasta(fasta_env), add = TRUE)

  expect_warning(
    result <- get_repetition_seq_info(
      chr = "chr_empty", pos = 1, ref = "G", alt = "G", # substring("G", 2) is ""
      fasta = fasta_env$fa_obj, mutation_type = "insertion"
    ),
    # First warning is for same length, if mutation_type logic changes, this might adapt.
    regexp = "Derived indel sequence is empty"
  )
  expect_null(result)

  # A more direct test for the empty indel sequence block:
  # Assume ref = "N", alt = "A", mutation_type = "insertion"
  # indel_seq becomes substring("A", 2) which is ""
  expect_warning(
    result_direct_empty <- get_repetition_seq_info(
      chr = "chr_empty", pos = 1, ref = "N", alt = "A", # alt="A" -> indel_seq = ""
      fasta = fasta_env$fa_obj, mutation_type = "insertion"
    ),
    regexp = "Derived indel sequence is empty"
  )
  expect_null(result_direct_empty)
})

test_that("Empty indel_seq from deletion (ref allele too short) returns NULL with warning", {
  sequences <- c(chr_empty_del = "GATTACA")
  fasta_env <- setup_test_fasta(sequences)
  on.exit(cleanup_test_fasta(fasta_env), add = TRUE)

  # ref="A", so substring(ref,2) is empty.
  expect_warning(
    result <- get_repetition_seq_info(
      chr = "chr_empty_del", pos = 1, ref = "A", alt = "N", # ref="A" -> indel_seq = ""
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
      chr = "chr_nonexistent", # This chromosome is not in fasta_env
      pos = 10,
      ref = "C",
      alt = "CAT",
      fasta = fasta_env$fa_obj,
      mutation_type = "insertion"
    ),
    # The regexp should match the warning generated by the modified function.
    regexp = "Error trying to retrieve sequence from FASTA at chr_nonexistent:11[.] Details: .*not found"
  )
  expect_null(result)
})

test_that("Function handles indel at the very start of sequence correctly", {
  sequences <- c(chr_start = "ATATAG") # No anchor base in this test sequence if pos=0
  # Let's assume pos=0 implies "before first base" conceptually
  # And genome is just the repeat part.
  fasta_env <- setup_test_fasta(sequences)
  on.exit(cleanup_test_fasta(fasta_env), add = TRUE)

  # If pos = 0, current_genome_pos = 1. Indel "AT"
  # Genome: A T A T A G
  # Pos:    1 2 3 4 5 6
  # Indel:  A T A T A (mismatch at G)
  result <- get_repetition_seq_info(
    chr = "chr_start", pos = 0, ref = "", alt = "AT", # Indel "AT", ref is empty anchor
    fasta = fasta_env$fa_obj, mutation_type = "insertion"
  )
  # Genome: A (pos 1) vs Indel T -> mismatch. Expected: 1
  expect_equal(result, 1)

  # Let's use a standard anchor: ref="N", alt="NAT" (indel "AT")
  # FASTA: N ATATATAG (where N is an actual base at pos 1)
  sequences2 <- c(chr_start_anchor = "NATATATAG")
  fasta_env2 <- setup_test_fasta(sequences2)
  on.exit(cleanup_test_fasta(fasta_env2), add = TRUE)

  result2 <- get_repetition_seq_info(
    chr = "chr_start_anchor", pos = 1, ref = "N", alt = "NAT",
    fasta = fasta_env2$fa_obj, mutation_type = "insertion"
  )
  expect_equal(result2, 9)
})
