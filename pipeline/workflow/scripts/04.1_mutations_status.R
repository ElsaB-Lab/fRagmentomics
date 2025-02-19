# @created: 02 Dec 24
# @modified: 02 Dec 24
# @author: Juliette Samaniego
#
# This script reads a TSV file, determines mutation status (MUT vs WT), and writes out the results.
# The mutation status logic is as follows:
# 1. If both reads cover the position of interest:
#    - For a SNV-type mutation: Only if both reads have the ALT base do we consider the fragment mutated (MUT). Otherwise, it is WT.
#    - For an INDEL-type mutation: Only if both reads show the indel do we consider the fragment mutated (MUT). Otherwise, it is WT.
# 2. If only one read covers the position:
#    - For a SNV-type mutation: If at least one read has the ALT base, we consider it MUT. Otherwise, WT.
#    - For an INDEL-type mutation: If at least one read shows the indel, we consider it MUT. Otherwise, WT.

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))

# Parse command line arguments =====================================================================
args <- commandArgs(trailingOnly = TRUE)

# Stop if not enough arguments
if (length(args) < 3) {
  stop("Requires 3 arguments: input_file, position, output_file.")
}

input_file    <- args[1]
pos           <- as.numeric(args[2])  # Convert position to numeric
output_file   <- args[3]

# Function definition ==============================================================================
DetermineMutationStatus <- function(inputfile, position, outputfile) {
  # Read the input TSV
  df <- read_tsv(
    file = inputfile,
    show_col_types = FALSE,
    progress = FALSE,
    col_types = cols(
      Base_Read_1    = col_character(),
      Base_Read_2    = col_character(),
      Alt            = col_character(),
      mutation_type  = col_character(),
      Startpos_read_1 = col_double(),
      Startpos_read_2 = col_double(),
      read_length1   = col_double(),
      read_length2   = col_double(),
      Indel_1        = col_double(),
      Indel_2        = col_double()
    ),
    quote = ""
  )

  # Variables
  Alt_value <- df$Alt[1]
  mutation_type <- df$mutation_type[1]

  # Prepare string patterns for SNVs (mutation_type == "mutation")
  # These patterns are used to detect which reads contain the ALT base.
  Base_read_mut_mut <- paste0(Alt_value, "_", Alt_value)      # read1=ALT, read2=ALT
  Base_read_mut_NA  <- paste0(Alt_value, "_NA")         # read1=ALT, read2=NA
  Base_read_NA_mut  <- paste0("NA_", Alt_value)         # read1=NA,  read2=ALT

  # Create a combined base for easy checking
  df <- df %>%
    mutate(Base_Read = paste(Base_Read_1, Base_Read_2, sep = "_"))

  # Compute coverage booleans for each read
  # TRUE if the read extends beyond the position of interest
  df <- df %>%
    mutate(
      covers_read1 = (Startpos_read_1 + Read_length_1) > position,
      covers_read2 = (Startpos_read_2 + Read_length_2) > position
    )

  # Create the Mutation_Status column using row-wise logic
  df$Mutation_Status <- sapply(seq_len(nrow(df)), function(i) {
    row_data <- df[i, ]
    both_cover <- row_data$covers_read1 & row_data$covers_read2

    if (mutation_type == "mutation") {
      # SNV logic
      if (both_cover) {
        # Both reads cover the position -> MUT only if both are ALT
        if (row_data$Base_Read == Base_read_mut_mut) {
          return("MUT")
        } else {
          return("WT")
        }
      } else {
        # One read or none covers the position -> MUT if at least one is ALT
        if (row_data$Base_Read %in% c(Base_read_mut_mut, Base_read_mut_NA, Base_read_NA_mut)) {
          return("MUT")
        } else {
          return("WT")
        }
      }
    } else {
      # Indel logic
      # We assume Indel_1 and Indel_2 are columns with 1 indicating the presence of the indel
      if (both_cover) {
        # Both reads cover -> MUT only if both reads have Indel
        if (row_data$Indel_1 == 1 && row_data$Indel_2 == 1) {
          return("MUT")
        } else {
          return("WT")
        }
      } else {
        # One read or none covers -> MUT if at least one read has Indel
        if (row_data$Indel_1 == 1 || row_data$Indel_2 == 1) {
          return("MUT")
        } else {
          return("WT")
        }
      }
    }
  })

  # Write the final output
  write_tsv(df, file = outputfile)
}

# Run the function =================================================================================
DetermineMutationStatus(
  inputfile   = input_file,
  position    = pos,
  outputfile  = output_file
)
