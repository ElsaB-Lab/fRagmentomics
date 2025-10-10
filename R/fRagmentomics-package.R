#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

#' fRagmentomics: Per-Fragment Analysis of cfDNA characteristics
#'
#' @description
#' A user-friendly R package that enables the characterization of each cfDNA
#' fragment overlapping one or multiple mutations of interest, starting from a
#' sequencing file containing aligned reads (BAM file). fRagmentomics supports
#' multiple mutation input formats (e.g., VCF, TSV, or string
#' “chr:pos:ref:alt” representation), accommodates one-based and zero-based
#' genomic conventions, handles mutation representation ambiguities, and accepts
#' any reference file and species in FASTA format.
#'
#' The core functionality is delivered through the `run_fRagmentomics()`
#' function, which implements a robust pipeline that:
#' \itemize{
#'   \item Accepts variant inputs in multiple formats (VCF, TSV, or string).
#'   \item Performs validation and normalization of variants into a canonical,
#'   left-aligned representation.
#'   \item Efficiently queries BAM files to retrieve all relevant fragments for
#'   each variant.
#'   \item Executes the analysis in parallel to ensure high performance on
#'   large datasets.
#'   \item For each fragment, it extracts a rich set of features, including its
#'   inferred size, mutational status (mutant, wild-type, or ambiguous),
#'   alignment details, and the sequences of its 5' and 3' ends.
#' }
#'
#' @seealso
#' The primary function for analysis: \code{\link{run_fRagmentomics}}
#'
#' @examples
#' # The main entry point of the package is the run_fRagmentomics() function.
#' ?run_fRagmentomics
#'
#' @docType package
#' @name fRagmentomics