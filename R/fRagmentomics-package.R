#' @title fRagmentomics: Characterization of cfDNA Fragments
#'
#' @description A user-friendly R package that enables the characterization
#'   of each cfDNA fragment overlapping one or multiple mutations of interest,
#'   starting from a sequencing file containing aligned reads (BAM file).
#'   fRagmentomics supports multiple mutation input formats (e.g., VCF, TSV,
#'   or string “chr:pos:ref:alt” representation), accommodates one-based and
#'   zero-based genomic conventions, handles mutation representation ambiguities,
#'   and accepts any reference file and species in FASTA format. For each cfDNA
#'   fragment, fRagmentomics outputs its size, its 3’ and 5’ sequences, and its
#'   mutational status.
#'
#' @author Killian Maudet <killian.maudet@gustaveroussy.fr> [aut, cre]
#' @author Yoann Pradat <yoann.pradat@gustaveroussy.fr> [aut, cre]
#' @author Juliette Samaniego <juliette.samaniego@gustaveroussy.fr> [aut]
#' @author Elsa Bernard <elsa.bernard@gustaveroussy.fr> [aut]
#'
#' @examples
#' # Example usage of fRagmentomics
#' # Load the package
#' library(fRagmentomics)
#'
#' # Example function call
#' # result <- analyze_fragments(bam_file = "path/to/your/file.bam",
#' #                             mutation_file = "path/to/your/mutations.vcf")
#'
#' @name fRagmentomics
NULL
