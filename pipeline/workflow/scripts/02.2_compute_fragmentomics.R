# @created: 08 Apr 24
# @modified: 08 Apr 24
# @authors: Yoann Pradat
#
# Run fRagmentomics R package

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(fRagmentomics))
suppressPackageStartupMessages(library(readr))

# functions ============================================================================================================

save_table <- function(table, filepath) {
  # save
  if (grepl(".gz$", filepath)) {
    write.table(table, gsub(".gz$", "", filepath), sep = "\t", row.names = F, quote = F)
    system(paste("gzip", gsub(".gz$", "", filepath)))
  } else {
    write.table(table, filepath, sep = "\t", row.names = F, quote = F)
  }
  cat("-table saved at", filepath, "\n")
}


main <- function(args) {
  sam_head <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual")
  df_sam <- read_tsv(args$sam, col_names = sam_head, show_col_types = F, progress = F)
  df_info <- build_fragments_info_table(
    df_sam = df_sam, sample_id = args$sample_id, chr = args$chr, pos = args$pos,
    ref = args$ref, alt = args$alt, del_info = args$del_info,
    mutation_type = args$mutation_type, n_cores = args$n_cores
  )
  save_table(df_info, args$out)
}


# run ==================================================================================================================

if (getOption("run.main", default = TRUE)) {
  parser <- ArgumentParser(description = "Run fRagmentomics package.")
  parser$add_argument("--sam", type = "character", help = "Path to input SAM file.")
  parser$add_argument("--sample_id", type = "character", help = "Sample_id of the mutation of interest.")
  parser$add_argument("--chr", type = "character", help = "Chromosome of the mutation of interest.")
  parser$add_argument("--pos", type = "integer", help = "Position of the mutation of interest.")
  parser$add_argument("--ref", type = "character", help = "Mutation reference allele.")
  parser$add_argument("--alt", type = "character", help = "Mutation alternatif allele.")
  parser$add_argument("--del_info", type = "character", help = "character vector sous la forme $pos_final,$rep_del")
  parser$add_argument("--mutation_type", type = "character", help = "mutation type in mutation, insertion, deletion")
  parser$add_argument("--n_cores", type = "integer", help = "Number of cores to be used.", default = 1)
  parser$add_argument("--out", type = "character", help = "Path to output table.")
  parser$add_argument("--log", type = "character", help = "Path to log file.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open = "wt")
  sink(log)
  sink(log, type = "message")

  print(args)
  main(args)
}
