pkgname <- "fRagmentomics"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "fRagmentomics-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('fRagmentomics')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("analyze_fragments")
### * analyze_fragments

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: analyze_fragments
### Title: Analyze fragments
### Aliases: analyze_fragments

### ** Examples

## Not run: 
##D # Load the package
##D library(fRagmentomics)
##D 
##D # Assuming you have your input files:
##D mut_file <- "path/to/your/mutation.tsv"
##D bam_file <- "path/to/your/alignment.bam"
##D fasta_file <- "path/to/your/reference.fasta"
##D output_path <- "path/to/your/results.tsv"
##D 
##D # Run the analysis on 4 cores
##D results_df <- analyze_fragments(
##D   mut = mut_file,
##D   bam = bam_file,
##D   fasta = fasta_file,
##D   output_file = output_path,
##D   n_cores = 4
##D )
##D 
##D # View the first few results
##D head(results_df)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("analyze_fragments", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fRagmentomics")
### * fRagmentomics

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fRagmentomics
### Title: fRagmentomics: Per-Fragment Analysis of cfDNA characteristics
### Aliases: fRagmentomics

### ** Examples

# The main entry point of the package is the analyze_fragments() function.

?analyze_fragments




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fRagmentomics", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
