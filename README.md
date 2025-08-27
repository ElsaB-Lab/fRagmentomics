# fRagmentomics <img src="man/figure/logo.svg" align="right" />

![Version](https://img.shields.io/badge/version-0.2.5-blue)
[![codecov](https://codecov.io/gh/ElsaB-Lab/fRagmentomics/graph/badge.svg?token=OMTSCRO7LJ)](https://codecov.io/gh/ElsaB-Lab/fRagmentomics)

## Overview

Plasma circulating cell-free DNA (cfDNA) analysis has transformed cancer care. The majority of cfDNA originates from hematopoietic cells, complicating the interpretation of ctDNA results in clinical practice, in the absence of matched white blood cells sequencing. Recent work has demonstrated that ctDNA fragments have distinct size distribution profiles and 5’/3’ end sequences compared to healthy cfDNA fragments.

Cependant, à l'heure actuelle, aucun outil publié ne standardsie une analyse approfondit des fragmentomics features associées avec le status mutationnel en prenant en compte les challenges spécifiques associés
aux cfDNA. En effet, pour n'en citer que quelqu'uns, il esxitse plusieurs approches du calcul de la fragment size lié. Simplement se basé sur les positions génomiques de début et de fin du fragment ne permet pas
d'avoir la taille effective car on ne prend pas en comtpe les indels au sein du fragment qui peuvent changer la taille du fragment. Il existe aussi plsuieurs représentations des mutations en fonction des callers ainsi
que des positions pour les indels.

**fRagmentomics** is an R package designed to be a standardised fragmeork which integrate cell-free DNA (cfDNA) fragment features (size, end sequences) with mutational status (tumor drivers or clonal hematopoiesis mutations) to support the interpretation of liquid biopsies.

Feature selection is based on the study by Wong *et al.* (2024):
**"Cell-free DNA from germline TP53 mutation carriers reflect cancer-like fragmentation patterns"**, *Nature Communications*.

---

## Installation

You can install the development version of **fRagmentomics** by activating the package environment from the root of the repository:

```bash
conda env create -f fRagmentomics.yaml
conda activate fRagmentomics
```

Then, open an R session and install the package with:

```r
devtools::install_github("ElsaB-Lab/fRagmentomics")
```

---

## Example Usage

Here is a basic example of how to call the main function:

```r
library(fRagmentomics)

process_fragmentomics(
  mut = "/path/to/mutation_file.tsv",
  bam = "/path/to/bam_file.bam",
  fasta = "/path/to/fasta_file.fa",
  sample_id = "Sample_Id",
  neg_offset_mate_search = -500,
  pos_offset_mate_search = 500,
  one_based = TRUE,
  flag_keep = 0x03,
  flag_remove = 0x900,
  report_tlen = TRUE,
  report_softclip = TRUE,
  report_5p_3p_bases_fragment = 5,
  tmp_folder = tempdir(),
  output_file = "/path/to/output_file.tsv",
  n_cores = 3
)
```

### Arguments

- **`mut`**: Mutation input. Accepts a proper VCF file, a TSV file with `CHROM`, `POS`, `REF`, `ALT` columns, or a string like `"chr:pos:ref:alt"`.
- **`bam`**: BAM file containing the fragment information for your sample.
- **`fasta`**: Reference FASTA file (must include a `.fai` index in the same directory).
- **`sample_id`**: Identifier for the sample.
- **`neg_offset_mate_search`**: Integer. Number of nucleotides to extend *upstream* from the position of interest when querying the BAM file.
- **`pos_offset_mate_search`**: Integer. Number of nucleotides to extend *downstream* from the position of interest.
- **`one_based`**: Logical. Set to `TRUE` if the BAM file uses one-based indexing.
- **`flag_keep`**: SAM flag for reads to keep. See [Picard Explain Flags](https://broadinstitute.github.io/picard/explain-flags.html).
- **`flag_remove`**: SAM flag for reads to exclude. See [Picard Explain Flags](https://broadinstitute.github.io/picard/explain-flags.html).
- **`report_tlen`**: Logical. Whether to report the TLEN (template length) of the fragment.
- **`report_softclip`**: Logical. Whether to report soft-clipped bases at the 5' and 3' ends of fragments.
- **`report_5p_3p_bases_fragment`**: Integer. Number of bases to report from both ends (5' and 3') of each fragment.
- **`tmp_folder`**: Path to a temporary directory for intermediate files.
- **`output_file`**: Path to the output TSV file.
- **`n_cores`**: Number of CPU cores to use for parallel processing.
