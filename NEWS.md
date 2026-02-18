# fRagmentomics NEWS

## fRagmentomics 0.99.9 (2026-02-18)

### OTHER CHANGES

- Preparing pre-release version.
- **Test**: Add test to the package.

## fRagmentomics 0.99.4 (2026-02-17)

### NEW FEATURES

- **Example Dataset**: Added a benchmark case study featuring a 15-bp in-frame deletion in *EGFR* exon 19 (lung cancer cfDNA) to demonstrate indel-specific fragment sizing and genotyping.
- **Reproducibility Scripts**: Included data generation and anonymization scripts for the test dataset in `inst/scripts/`.
- **Developer Tools**: Added a pre-commit configuration to standardize the development workflow.

### OTHER CHANGES

- **Documentation**:
  - Updated README and vignette to incorporate the new case study and illustrative figures.
- **Testing**: Updated the test suite to include coverage for the new *EGFR* case study.

## fRagmentomics 0.99.3 (2026-02-16)

### OTHER CHANGES

- **QC Logic**: Improved `process_fragment_reads_QC` with updated default parameters and enhanced messaging for failed fragments.
- **Code Refactoring**: Factorized internal code and replaced `stop()` calls in `extract_fragment_features` to improve robustness.

## fRagmentomics 0.99.2 (2025-10-28)

### OTHER CHANGES

- Improved `future.apply::future_lapply` parallelization by setting properly the `chunk` parameter and reduced exported objects per worker.

## fRagmentomics 0.99.1 (2025-10-22)

### NEW FEATURES

- Added `flag_bam_list` instead of `TLEN` (`default: FALSE`). Parameter to include all BAM fields in the output.
- Added `Position_3p` to the output. `Position_3p` represents the last aligned position of the 3′ read of the fragment.
- Replaced rbind in the main function with `data.table::rbindlist` to make it more consistent and faster.

### OTHER CHANGES

- Change `future.apply::future_lapply` parameter to remove a warning.

---

## fRagmentomics 0.99.0 (2025-10-09)

### NEW FEATURES

- Added `verbose` (default `FALSE`) parameter to remove messages.

## fRagmentomics 0.2.9 (2025-10-07)

### NEW FEATURES

- Change fragment if DIS with potentially compatible with DIS.

### OTHER CHANGES

- Update vignette.
- Fix some display errors on graph functions.

---

## fRagmentomics 0.2.8 (2025-08-27)

### NEW FEATURES

- Added new plotting functions for in-depth fragment analysis:
  - `plot_size_distribution()` to visualize fragment size distributions.
  - `plot_freq_barplot()` to display nucleotide frequencies in end motifs.
  - `plot_ggseqlogo_meme()` for generating sequence logo plots of end motifs using `ggseqlogo`.
  - `plot_motif_barplot()` for generating end motifs in barplot.

### OTHER CHANGES

- The package has been updated to meet all Bioconductor submission requirements.

---

## fRagmentomics 0.2.7 (2025-08-15)

### IMPROVEMENTS

- Improved fragment size calculation to more accurately account for insertions and deletions within the overlapping section of read pairs.
- Added the `remove_softclip` parameter to trim soft-clipped bases from the ends of fragments.

---

## fRagmentomics 0.2.6 (2025-07-30)

### IMPROVEMENTS

- Upgraded BAM flag selection for more precise read filtering.
- Fragments with reads aligned in the same orientation are now automatically removed.
- The definition of 5' and 3' reads is now based on strand information rather than genomic position for greater accuracy.

---

## fRagmentomics 0.2.5 (2025-07-10)

### NEW FEATURES

- Introduced parallel processing (`future`) and a progress bar (`progressr`) to significantly speed up analysis on multi-core systems.

### IMPROVEMENTS

- Improved mutation classification for reads with soft-clipping at their ends, preventing potential false negatives. Cases previously classified as `WT` are now correctly labeled `AMB`.
- Added an `Input_Mutation` column to the output to retain the original mutation information before normalization.

---

## fRagmentomics 0.2.4 (2025-06-25)

### BUG FIXES

- Fixed a critical bug related to incorrect fetching of the reference sequence.
- Improved handling of ambiguous cases and multi-nucleotide variants (MNVs).
- Optimized the SNV detection algorithm.

---

## fRagmentomics 0.2.3 (2025-06-12)

### IMPROVEMENTS

- Improved performance by replacing per-fragment FASTA requests with a single, larger request per mutation, resulting in a **1.5x to 2x speed increase**.

---

## fRagmentomics 0.2.2 (2025-05-20)

### IMPROVEMENTS

- Introduced more descriptive mutation statuses (e.g., `"MUT but potentially larger MNV"`, `"OTH (DEL)"`).
- Updated the reporting of `BASE_5p` and `BASE_3p` to be more comprehensive around the variant position.
- Finalized the definitions for `Fragment_Status_Simple` and `Fragment_Status_Detail` columns.

### BUG FIXES

- Corrected mutation status assignment for reads with hard or soft-clipping near the variant.
- The `get_index_aligning_with_pos()` function now correctly handles deletions at the position of interest.

### OTHER CHANGES

- Renamed parameter `cigar_free_mode` to `cigar_free_indel_match` for clarity.

---

## fRagmentomics 0.2.1 (2025-05-01)

### BUG FIXES

- Corrected the Variant Allele Frequency (VAF) calculation.

### OTHER CHANGES

- Harmonized all output column names to use "Upper_Snake_Case" for consistency.

---

## fRagmentomics 0.2.0 (2025-04-15)

### NEW FEATURES

- Implemented a new, alternative algorithm for determining mutation status that is independent of the CIGAR string. This significantly improves INDEL detection and is controllable via the `cigar_free_indel_match` parameter.

### OTHER CHANGES

- Major code refactoring for clearer function names and better maintainability.

---

## fRagmentomics 0.1.0 (2025-03-01)

### NEW FEATURES

- Initial beta version of the package.
- Core functionality for extracting key fragmentomic features, including **end motifs** and **fragment size**.
- Initial implementation of fragment mutational status classification (e.g., `WT`, `MUT`, `AMB`, and more complex states)
