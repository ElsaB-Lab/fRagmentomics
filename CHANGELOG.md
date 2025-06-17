# fRagmentomics 0.2.2

- Correct mutation read status to handle incidental mutations arising from hard/soft-clipping.
- Correct get_index_aligning_with_pos to return -2 when the at the position of interest the read bears a deletion.
- Implement more verbose mutation statuses
    - "MUT but potentially larger MNV"
    - "OTH (DEL)" if the read bears a deletion at the position of interest and we are looking for a SNV/MNV
    - "OTH (DEL) for indel case - Need review" if the read bears a deletion at the position of interest for an INDE:.
- Rename `cigar_free_mode` to `cigar_free_indel`. This option, when activated, serves to search for an INDEL if it is
  not found by the function `search_for_indel_in_cigar`.
- Update how `BASE_5p` and `BASE_3p` are reported.  The idea is to report the bases aligning with pos, [inserted bases],
  pos+1, [inserted bases], ..., pos+nchar(alt) and stop reporting bases once we reach nchar(alt).
- Implement final Fragment_Status_Simple and Fragment_Status_Detail values.
- No graphic representation
- No pipeline integration.

# fRagmentomics 0.2.1

- Correct VAF calculation
- Harmonize column names to first-letter upper case and underscores.
- No graphic representation
- No pipeline integration.

# fRagmentomics 0.2.0

- Implement an alternative algorithm for determining the mutation status of a fragment. This algorithm relies on the
  idea that by comparing the read sequence to the wild-type reference and the "mutated" reference sequences we can
  determine the mutation status for any mutation type (SNV, MNV, INS, DEL) and regardless of the mutation
  representation.
  The downside of this algorithm is its potential for false positives. To mitigate this, a cigar_free_mode option has
  been added to allow for two running modes. This option only influences the mutation status for INDELS.
    1. `cigar_free_mode=FALSE` then the mutation status relies on the code from `0.1.0` to determine the mutation status
       from the CIGAR. If the mutation is not found in the CIGAR but comparison to wild-type/mutated ref sequences
       detects the mutation, then the mutation status is to specific labels to track these cases.
    1. `cigar_free_mode=TRUE` then the mutation status is determined completely independently of the CIGAR.
- The code was refactored to have clearer function names wherever possible.
- No graphic representation
- No pipeline integration.

# fRagmentomics 0.1.0

- Beta version.
- Extraction of fragmentomic features: **end motifs, fragment size**.
- Extract fragment status:
  - **WT**: Both reads WT or only one if only one covers the position of interest,
  - **MUT**: Both reads MUT or only one if only one covers the position of interest,
  - **WT_but_other_read_MUT**: One read MUT and the other one WT (both cover),
  - **WT_but_other_read_mut_with_other_alt**: One read WT and the other one MUT with an other mutation (both cover),
  - **MUT_but_other_read_mut_with_other_alt**: One read MUT and the other one MUT with an other mutation (both cover),
  - **Other_MUT**: One read is covering and it's not ref or alt,
  - **Error_both_read_mut_with_other_alt**: Both reads cover the position but both are not ref or alt. Error because the mutation is supposed to be found in the sample.
- For indels, if the position (nucleotide before the indel) is the last nucleotide of the read, the read will be considered as covering.
- No graphic representation
- No pipeline integration.
