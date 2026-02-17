#' Compare a read sequence against reference and alternate alleles
#'
#' @description A core comparison utility that classifies a read sequence by
#' checking its compatibility with the expected wild-type (WT) and mutant (MUT)
#' sequences.
#'
#' @param read_seq Base sequence of the read
#' @param ref_seq_wt Base sequence of the wild-type ref
#' @param ref_seq_mut Base sequence of the mutated ref
#' @param compare_len_wt Number of bases to be compared in read vs wild-type ref
#'   comparison
#' @param compare_len_mut Number of bases to be compared in read vs
#'   mutated ref comparison
#' @return  a character describing the mutation status of the read
#'
#' @keywords internal
compare_read_to_ref_wt_and_mut <- function(
  read_seq, ref_seq_wt, ref_seq_mut,
  compare_len_wt, compare_len_mut
) {
  # build the sequences to be compared
  ref_seq_wt_sub <- substr(ref_seq_wt, 1, compare_len_wt)
  ref_seq_mut_sub <- substr(ref_seq_mut, 1, compare_len_mut)
  read_seq_wt_sub <- substr(read_seq, 1, compare_len_wt)
  read_seq_mut_sub <- substr(read_seq, 1, compare_len_mut)

  # check if read compatible with wild-type ref
  match_ref_wt <- ref_seq_wt_sub == read_seq_wt_sub
  # check if read compatible with mutated ref
  match_ref_mut <- ref_seq_mut_sub == read_seq_mut_sub

  if (match_ref_wt && match_ref_mut) {
    return("AMB")
  } else if (match_ref_wt) {
    return("WT")
  } else if (match_ref_mut) {
    return("MUT")
  } else {
    return("OTH")
  }
}
