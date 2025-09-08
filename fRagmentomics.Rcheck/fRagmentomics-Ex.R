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

# --- 1. Locate Example Files ---
# The package includes small example files to demonstrate its functionality.
# We locate them using system.file().
mut_file <- system.file(
  "extdata", "mutations_cfdna-test-01_chr1_27433000_27435000.tsv",
  package = "fRagmentomics"
)
bam_file <- system.file(
  "extdata", "cfdna-test-01_chr1_27433000_27435000.bam",
  package = "fRagmentomics"
)
fasta_file <- system.file(
  "extdata", "hg19_chr1_27433000_27435000.fa",
  package = "fRagmentomics"
)

# --- 2. Run the Analysis ---
# This single call runs the full analysis pipeline on the example data.
# The output file is written to a temporary location to avoid cluttering
# the working directory. We use n_cores = 1L for examples.
results <- analyze_fragments(
  mut = mut_file,
  bam = bam_file,
  fasta = fasta_file,
  sample_id = "cfdna-test-01",
  output_folder = tempdir(),
  n_cores = 1L
)

# --- 3. View the Results ---
# Print the first few rows of the output data frame to see the results.
print(head(results))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("analyze_fragments", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_freq_barplot")
### * plot_freq_barplot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_freq_barplot
### Title: Plot Overall Nucleotide Frequency
### Aliases: plot_freq_barplot

### ** Examples

## --- Create a dataset for demonstration ---
# Set a seed for reproducibility
set.seed(42)

# Helper function to generate random DNA sequences with a bias
generate_biased_dna <- function(n_seq, len, prob) {
    bases <- c("A", "C", "G", "T")
    replicate(n_seq, paste(sample(bases, len, replace = TRUE, prob = prob), collapse = ""))
}

# Create 50 "MUT" fragments with a high proportion of 'C'
df_mut <- data.frame(
    Fragment_Bases_5p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
    Fragment_Bases_3p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
    Fragment_Status_Simple = "MUT"
)

# Create 50 "WT" fragments with a high proportion of 'G'
df_wt <- data.frame(
    Fragment_Bases_5p = generate_biased_dna(50, 10, prob = c(0.15, 0.15, 0.5, 0.2)),
    Fragment_Bases_3p = generate_biased_dna(50, 10, prob = c(0.15, 0.15, 0.5, 0.2)),
    Fragment_Status_Simple = "WT"
)

# Combine into a single dataframe
example_df <- rbind(df_mut, df_wt)

## --- Function Calls ---

# 1. Default plot: Compares MUT vs. WT groups for 3-mers
#    from both the 5' and 3' ends.
p1 <- plot_freq_barplot(example_df)
print(p1)

# 2. Customized plot: Analyzes only the first nucleotide ('motif_size = 1')
#    of the 5' end ('motif_type = "Start"') using custom colors.
p2 <- plot_freq_barplot(
    df_fragments = example_df,
    motif_type = "Start",
    motif_size = 1,
    colors_z = c("MUT" = "#d95f02", "WT" = "#1b9e77")
)
print(p2)

# 3. Ungrouped plot: Analyzes the overall nucleotide frequency
#    across all fragments combined.
p3 <- plot_freq_barplot(example_df, col_z = NULL)
print(p3)

# 4. Plot with a subset of groups: If you had more than two groups
#    (e.g., "MUT", "WT", "AMB"), you could select specific ones to plot.
p4 <- plot_freq_barplot(
    df_fragments = example_df,
    vals_z = c("MUT", "WT")
)
print(p4)

# 5. Save the default plot to a temporary folder.
# plot_freq_barplot(
#   df_fragments = example_df,
#   sample_id = "test01_freq",
#   output_folder = tempdir()
# )

# 6. Save a customized plot with specific dimensions.
# plot_freq_barplot(
#   df_fragments = example_df,
#   motif_type = "Start",
#   motif_size = 1,
#   sample_id = "test02_freq_custom",
#   output_folder = tempdir(),
#   ggsave_params = list(width = 7, height = 5, units = "in")
# )




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_freq_barplot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_motif_barplot")
### * plot_motif_barplot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_motif_barplot
### Title: Plot 3-base motif proportions with various representations
### Aliases: plot_motif_barplot

### ** Examples

## --- Create a dataset for demonstration ---
# Set a seed for reproducibility
set.seed(42)

# Helper function to generate random DNA sequences with a bias
generate_biased_dna <- function(n_seq, len, prob) {
    bases <- c("A", "C", "G", "T")
    replicate(n_seq, paste(sample(bases, len, replace = TRUE, prob = prob), collapse = ""))
}

# Create 50 "MUT" fragments with a high proportion of motifs starting with 'C'
df_mut <- data.frame(
    Fragment_Bases_5p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
    Fragment_Bases_3p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
    Fragment_Status_Simple = "MUT"
)

# Create 50 "WT" fragments with a high proportion of motifs starting with 'G'
df_wt <- data.frame(
    Fragment_Bases_5p = generate_biased_dna(50, 10, prob = c(0.15, 0.15, 0.5, 0.2)),
    Fragment_Bases_3p = generate_biased_dna(50, 10, prob = c(0.15, 0.15, 0.5, 0.2)),
    Fragment_Status_Simple = "WT"
)

# Combine into a single dataframe
example_df <- rbind(df_mut, df_wt)

## --- Function Calls for Each Representation ---

# 1. Hierarchical Plot (representation = "split_by_base")
# This is the default. It creates nested facets for each base position.
p1 <- plot_motif_barplot(
    df_fragments = example_df,
    representation = "split_by_base"
)
print(p1)

# You can also filter this plot to show only motifs starting with certain bases.
p1_filtered <- plot_motif_barplot(
    df_fragments = example_df,
    representation = "split_by_base",
    motif_start = c("C", "G")
)
print(p1_filtered)

# 2. Differential Plot (representation = "differential")
# This shows the log2 fold change in motif proportions between two groups.
# It requires exactly two groups specified in 'vals_z'.
p2 <- plot_motif_barplot(
    df_fragments = example_df,
    representation = "differential",
    vals_z = c("MUT", "WT")
)
print(p2)

# 3. Side-by-side Motif Plot (representation = "split_by_motif")
# This creates a more traditional bar plot with motifs on the x-axis and
# bars for each group shown side-by-side.
p3 <- plot_motif_barplot(
    df_fragments = example_df,
    representation = "split_by_motif"
)
print(p3)

# 4. Save the default hierarchical plot.
# plot_motif_barplot(
#   df_fragments = example_df,
#   sample_id = "test01_hierarchical",
#   output_folder = tempdir()
# )

# 5. Save the differential plot with custom dimensions.
# plot_motif_barplot(
#   df_fragments = example_df,
#   representation = "differential",
#   vals_z = c("MUT", "WT"),
#   sample_id = "test02_differential",
#   output_folder = tempdir(),
#   ggsave_params = list(width = 12, height = 8, units = "in")
# )




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_motif_barplot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_qqseqlogo_meme")
### * plot_qqseqlogo_meme

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_qqseqlogo_meme
### Title: Plot sequence motif composition
### Aliases: plot_qqseqlogo_meme

### ** Examples

## --- Create a dataset for demonstration ---
# Set a seed for reproducibility
set.seed(42)

# Helper function to generate random DNA sequences with a bias
generate_biased_dna <- function(n_seq, len, prob) {
    bases <- c("A", "C", "G", "T")
    replicate(n_seq, paste(sample(bases, len, replace = TRUE, prob = prob), collapse = ""))
}

# Create 50 "MUT" fragments with a high proportion of 'C' at the ends
df_mut <- data.frame(
    Fragment_Bases_5p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
    Fragment_Bases_3p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
    Fragment_Status_Simple = "MUT"
)

# Create 50 "WT" fragments with a high proportion of 'G' at the ends
df_wt <- data.frame(
    Fragment_Bases_5p = generate_biased_dna(50, 10, prob = c(0.15, 0.15, 0.5, 0.2)),
    Fragment_Bases_3p = generate_biased_dna(50, 10, prob = c(0.15, 0.15, 0.5, 0.2)),
    Fragment_Status_Simple = "WT"
)

# Combine into a single dataframe
example_df <- rbind(df_mut, df_wt)

## --- Function Calls ---

# 1. Default plot: Shows a 3-base motif from both 5' and 3' ends,
#    separated by a dash, for each group ("MUT" and "WT").
p1 <- plot_qqseqlogo_meme(example_df)
print(p1)

# 2. Analyze a longer, single-end motif: Shows a 5-base motif
#    from only the 5' end ('motif_type = "Start"').
p2 <- plot_qqseqlogo_meme(
    df_fragments = example_df,
    motif_type = "Start",
    motif_size = 5
)
print(p2)

# 3. Customizing colors: Use a named RColorBrewer palette.
#    Note the separator "-" is not a nucleotide and won't be colored.
p3 <- plot_qqseqlogo_meme(
    df_fragments = example_df,
    colors_z = "Set1"
)
print(p3)

# You can also provide a named vector for full control over colors.
custom_cols <- c("A" = "#1B9E77", "C" = "#D95F02", "G" = "#7570B3", "T" = "#E7298A")
p4 <- plot_qqseqlogo_meme(
    df_fragments = example_df,
    motif_type = "Start",
    colors_z = custom_cols
)
print(p4)

# 4. Ungrouped plot: Analyzes all fragments together as a single group.
p5 <- plot_qqseqlogo_meme(example_df, col_z = NULL)
print(p5)

# 5. Save plot with default settings.
# plot_qqseqlogo_meme(
#   df_fragments = example_df,
#   sample_id = "test01_motif",
#   output_folder = tempdir()
# )

# 6. Save plot with custom dimensions.
# plot_qqseqlogo_meme(
#   df_fragments = example_df,
#   sample_id = "test02_motif_custom",
#   output_folder = tempdir(),
#   ggsave_params = list(width = 15, height = 10, units = "cm")
# )




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_qqseqlogo_meme", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_size_distribution")
### * plot_size_distribution

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_size_distribution
### Title: Plot Fragment Size Distribution
### Aliases: plot_size_distribution

### ** Examples

## --- Create a dataset for demonstration ---
# Set a seed for reproducibility
set.seed(42)

# Generate fragment sizes for two groups with different distributions
# "MUT" group: N=100, shorter fragments
mut_sizes <- rnorm(100, mean = 150, sd = 20)

# "WT" group: N=150, centered around the mononucleosome peak
wt_sizes <- rnorm(150, mean = 170, sd = 25)

# Add some larger, dinucleosomal fragments to both groups
di_nuc_sizes <- rnorm(30, mean = 330, sd = 30)

# Combine into a single dataframe
example_df_size <- data.frame(
  Fragment_Size = c(mut_sizes, wt_sizes, di_nuc_sizes),
  Fragment_Status_Simple = c(
    rep("MUT", 100),
    rep("WT", 150),
    sample(c("MUT", "WT"), 30, replace = TRUE)
  )
)
# Ensure all fragment sizes are positive
example_df_size <- example_df_size[example_df_size$Fragment_Size > 0, ]

## --- Plotting Examples ---

# 1. Default plot: A grouped density plot with nucleosome peaks shown.
p1 <- plot_size_distribution(example_df_size)
print(p1)

# 2. Histogram plot: Show distributions as histograms instead of density curves.
#    We add transparency (alpha) so overlapping bars are visible.
p2 <- plot_size_distribution(
  df_fragments = example_df_size,
  show_histogram = TRUE,
  show_density = FALSE,
  histo_args = list(alpha = 0.6)
)
print(p2)

# 3. Combined plot: Overlay both density curves and histograms.
p3 <- plot_size_distribution(
  df_fragments = example_df_size,
  show_histogram = TRUE,
  show_density = TRUE,
  histo_args = list(alpha = 0.4)
)
print(p3)

# 4. Ungrouped and customized plot: Analyze all fragments together,
#    zoom in on the x-axis, and hide the nucleosome peak lines.
p4 <- plot_size_distribution(
  df_fragments = example_df_size,
  col_z = NULL,
  x_limits = c(50, 400),
  show_nuc_peaks = FALSE
)
print(p4)

# 5. Save plot with default settings (8x6 inches).
# plot_size_distribution(
#   df_fragments = example_df_size,
#   sample_id = "test01",
#   output_folder = tempdir()
# )

# 6. Save plot with custom dimensions in centimeters.
# plot_size_distribution(
#   df_fragments = example_df_size,
#   sample_id = "test02_custom",
#   output_folder = tempdir(),
#   ggsave_params = list(width = 20, height = 10, units = "cm", bg = "white")
# )




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_size_distribution", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
