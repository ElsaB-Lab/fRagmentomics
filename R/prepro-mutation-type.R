define_mutation_status <- function(ref, alt) {
  # Remove any underscores from ref and alt
  ref <- gsub("_", "", ref)
  alt <- gsub("_", "", alt)
  ref <- ifelse(ref == "NA", "", ref)
  alt <- ifelse(alt == "NA", "", alt)
  
  # Handle cases where ref or alt contains "NA"
  if (ref == "NA" && alt != "NA") {
    return("Insertion")
  } else if (alt == "NA" && ref != "NA") {
    return("Deletion")
  }
  
  # Determine mutation type
  if (nchar(ref) == nchar(alt)) {
    if (nchar(ref) == 1) {
      return("SNV")  # Single Nucleotide Variant
    } else {
      return("MNP")  # Multi-Nucleotide Polymorphism
    }
  } else if (nchar(ref) > nchar(alt)) {
    return("deletion")  # Deletion event
  } else {
    return("insertion")  # Insertion event
  }
}