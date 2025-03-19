suppressMessages(library(readr))

# read sam file
sam_head <- c("qname","flag","rname","pos","mapq","cigar","rnext","pnext","tlen","seq","qual")
sam_path <- system.file("testdata", "chr17_examples.sorted.sam", package="fRagmentomics")
df_sam <- readr::read_tsv(sam_path, col_names=sam_head, show_col_types=F, progress=F)

cat(paste("tempdir located at", tempdir(), "\n"))
