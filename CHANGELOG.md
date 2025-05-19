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
