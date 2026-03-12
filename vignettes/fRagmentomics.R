## ----setup, include=FALSE-----------------------------------------------------
use_ragg <- requireNamespace("ragg", quietly = TRUE)
knitr::opts_chunk$set(
  fig.width = 10,
  fig.height = 6,
  out.width = "100%",
  fig.align = "center",
  fig.retina = 2,
  dev = if (use_ragg) "ragg_png" else "png",
  dpi = 120
)

## ----echo=FALSE, results="hide", warning=FALSE--------------------------------
suppressPackageStartupMessages({
  library(fRagmentomics)
})

# Locate the example files bundled with the package
mut_file <- system.file(
  "extdata/mutation",
  "cfdna-egfr-del_chr7_55241864_55243064_10k.mutations.tsv",
  package = "fRagmentomics"
)
bam_file <- system.file(
  "extdata/bam",
  "cfdna-egfr-del_chr7_55241864_55243064_10k.bam",
  package = "fRagmentomics"
)
fasta_file <- system.file(
  "extdata/fasta",
  "hg19_chr7_55231864_55253064.fa",
  package = "fRagmentomics"
)

# Print the file paths to confirm they were found
mut_file
bam_file
fasta_file

# Detect whether bcftools is available on the system PATH
has_bcftools <- nzchar(Sys.which("bcftools"))
if (!has_bcftools) {
  message(
    "bcftools not found on PATH; normalization examples will run without it."
  )
}

## ----run_analysis_basic, message=TRUE, warning=TRUE---------------------------
if (!has_bcftools) {
  message("bcftools not detected; running without external normalization.")
}

# Run the full analysis pipeline with default settings
df_frag_basic <- run_fRagmentomics(
  mut = mut_file,
  bam = bam_file,
  fasta = fasta_file,
  sample_id = "cfdna-egfr-del",
  apply_bcftools_norm = has_bcftools,
  verbose = FALSE,
  n_cores = 1
)

# View the dimensions and the first few rows
cat("Dimensions of basic results:", dim(df_frag_basic), "\n")
head(df_frag_basic)

## ----run_analysis_advanced, message=TRUE, warning=TRUE------------------------
if (!has_bcftools) {
  message(
    "bcftools not detected; running advanced workflow without external normalization."
  )
}

# Run the analysis with custom parameters
df_frag_advanced <- run_fRagmentomics(
  mut = mut_file,
  bam = bam_file,
  fasta = fasta_file,
  sample_id = "cfdna-egfr-del",
  neg_offset_mate_search = -500,
  pos_offset_mate_search = 500,
  flag_bam_list = list(
    isProperPair = TRUE,
    isUnmappedQuery = FALSE,
    hasUnmappedMate = FALSE,
    isSecondaryAlignment = FALSE
  ),
  report_bam_info = TRUE,
  report_softclip = TRUE,
  retain_fail_qc = TRUE,
  report_5p_3p_bases_fragment = 3,
  apply_bcftools_norm = has_bcftools,
  verbose = FALSE,
  n_cores = 1
)

# Observe that stricter filtering may result in fewer fragments
cat("Dimensions of advanced results:", dim(df_frag_advanced), "\n")

# View the newly added columns
head(df_frag_advanced[, c(
  "Fragment_Id",
  "TLEN",
  "Nb_Fragment_Bases_Softclip_5p",
  "Fragment_Bases_5p"
)])

## ----dtaframe_output----------------------------------------------------------
# Dataframe output
df_frag_basic

## ----count_status_simple------------------------------------------------------
table(df_frag_advanced$Fragment_Status_Simple)

## ----count_status_detail------------------------------------------------------
table(df_frag_advanced$Fragment_Status_Detail)

## ----summary_size-------------------------------------------------------------
summary(df_frag_advanced$Fragment_Size)

## ----extract_vaf--------------------------------------------------------------
unique(df_frag_advanced$VAF)

