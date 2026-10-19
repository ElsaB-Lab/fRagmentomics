#' Analyze fragments overlapping genomic regions
#'
#' @description
#' End-to-end pipeline for extracting fragment-level features from a BED file
#' of genomic regions. Unlike `run_fRagmentomics_variant()`, this function does
#' not require a FASTA file and does not compute mutation status or VAF: it
#' reports the fragment features and the region each fragment overlaps.
#'
#' @details
#' For each region, the function queries the BAM file, retrieves complete
#' fragments with at least one read overlapping the region, and processes each
#' fragment in parallel. Results from all regions are aggregated into a single
#' data frame. Regions without any overlapping read are skipped with a warning.
#'
#' A heuristic warning lists all regions wider than 10000 bp. Wide regions are
#' not rejected: the threshold is only a signal that the output may aggregate
#' fragments with less specific coordinates.
#'
#' Fragment selection is based on overlap with the region: a fragment is
#' reported as soon as one of its reads shares at least one base with the
#' region. This convention is size-dependent, because for a fixed region width
#' longer fragments are more likely to overlap the region. The fragment size
#' distribution of a region is therefore enriched in long fragments compared
#' with the fragments originating from that locus, and the strength of this
#' enrichment depends on the region width. Size distributions should be
#' compared between regions of comparable width, or the selection criterion
#' should be changed (requiring a fragment end to fall inside the region is
#' size-neutral, while requiring full containment favours short fragments). A
#' fragment overlapping several regions is reported once per region.
#'
#' @inheritParams run_fRagmentomics_variant
#' @param region Path to a `.bed` or `.bed.gz` file. BED coordinates are
#'   0-based and half-open; they are converted internally to 1-based inclusive
#'   coordinates.
#'
#' @return A `DataFrame` (S4Vectors) containing extracted fragment-level
#'   information, or `NULL` when `output_path` is provided and the table is
#'   written to disk.
#'
#' @importFrom data.table rbindlist
#' @importFrom utils write.table
#' @importFrom future plan multisession sequential availableCores
#' @importFrom future.apply future_lapply
#' @importFrom S4Vectors DataFrame
#'
#' @export
#' @examples
#' # --- 1. Locate Example Files ---
#' region_file <- system.file(
#'   "extdata/bed", "exemple-regions.bed",
#'   package = "fRagmentomics"
#' )
#' bam_file <- system.file(
#'   "extdata/bam", "cfdna-egfr-del_chr7_55241864_55243064_10k.bam",
#'   package = "fRagmentomics"
#' )
#'
#' # --- 2. Run the Analysis ---
#' results <- run_fRagmentomics_region(
#'   region = region_file,
#'   bam = bam_file,
#'   sample_id = "cfdna-egfr-del",
#'   n_cores = 1L
#' )
#'
#' # --- 3. View the Results ---
#' head(results)
#'
run_fRagmentomics_region <- function(
  region, bam, sample_id = NA_character_, neg_offset_mate_search = -600,
  pos_offset_mate_search = 600, flag_bam_list = list(
    isPaired = TRUE,
    isProperPair = NA, isUnmappedQuery = FALSE, hasUnmappedMate = FALSE, isMinusStrand = NA,
    isMateMinusStrand = NA, isFirstMateRead = NA, isSecondMateRead = NA, isSecondaryAlignment = FALSE,
    isSupplementaryAlignment = FALSE, isNotPassingQualityControls = NA, isDuplicate = NA
  ),
  report_bam_info = FALSE, report_softclip = FALSE, report_5p_3p_bases_fragment = 5,
  remove_softclip = FALSE, retain_fail_qc = FALSE, output_path = NA_character_,
  verbose = FALSE, n_cores = 1
) {
  neg_offset_mate_search <- as.integer(neg_offset_mate_search)
  pos_offset_mate_search <- as.integer(pos_offset_mate_search)
  report_5p_3p_bases_fragment <- as.integer(report_5p_3p_bases_fragment)
  n_cores <- as.integer(n_cores)

  check_parameters_region(
    region, bam, sample_id, neg_offset_mate_search, pos_offset_mate_search,
    flag_bam_list, report_bam_info, report_softclip,
    report_5p_3p_bases_fragment, remove_softclip, retain_fail_qc,
    output_path, verbose, n_cores
  )

  df_regions <- read_bed(region)

  wide_regions <- df_regions$region_id[df_regions$width > 10000]
  if (length(wide_regions) > 0) {
    warning(sprintf(
      paste0(
        "The following region(s) are wider than 10000 bp: %s. ",
        "This threshold is heuristic; wide regions may aggregate less specific fragments."
      ),
      paste(wide_regions, collapse = ", ")
    ), call. = FALSE, immediate. = TRUE)
  }

  final_n_cores <- min(n_cores, future::availableCores())

  if (final_n_cores > 1L) {
    old_plan <- future::plan(future::multisession, workers = final_n_cores)
  } else {
    old_plan <- future::plan(future::sequential)
  }
  on.exit(future::plan(old_plan), add = TRUE)

  if (verbose) {
    message(sprintf("Processing %d region(s).", nrow(df_regions)))
  }

  results_list <- vector("list", nrow(df_regions))

  for (i in seq_len(nrow(df_regions))) {
    chr_region <- df_regions$chrom[i]
    start_region <- df_regions$start[i]
    end_region <- df_regions$end[i]
    region_id <- df_regions$region_id[i]

    if (verbose) {
      message(sprintf(
        "[%d/%d] Processing region %s (%s:%d-%d)",
        i, nrow(df_regions), region_id, chr_region, start_region, end_region
      ))
    }

    df_sam <- read_bam(
      bam, chr_region, start_region, end_region, neg_offset_mate_search,
      pos_offset_mate_search, flag_bam_list
    )

    if (is.null(df_sam)) {
      warning(sprintf(
        "No read covers the region of interest for the region %s. Skipping.",
        region_id
      ), immediate. = TRUE, call. = FALSE)
      next
    }

    fragment_names <- unique(df_sam$QNAME)

    number_chunks <- 50L
    chunk_size <- max(1L, ceiling(length(fragment_names) / number_chunks))

    chunk_indices <- ceiling(seq_along(fragment_names) / chunk_size)
    list_fragment_chunks <- split(fragment_names, chunk_indices)

    list_df_sam_chunks <- lapply(
      list_fragment_chunks,
      function(fragment_chunk) {
        df_sam[df_sam$QNAME %in% fragment_chunk, , drop = FALSE]
      }
    )

    list_results_all_chunks <- future.apply::future_lapply(
      seq_along(list_fragment_chunks),
      FUN = function(idx) {
        fragment_chunk <- list_fragment_chunks[[idx]]
        df_sam_chunk <- list_df_sam_chunks[[idx]]

        lapply(
          fragment_chunk,
          function(fragment_name) {
            extract_region_fragment_features(
              df_sam = df_sam_chunk,
              fragment_name = fragment_name,
              sample_id = sample_id,
              chr = chr_region,
              region_id = region_id,
              region_start = start_region,
              region_end = end_region,
              report_bam_info = report_bam_info,
              report_softclip = report_softclip,
              report_5p_3p_bases_fragment = report_5p_3p_bases_fragment,
              remove_softclip = remove_softclip
            )
          }
        )
      }
    )

    list_results_all_fragments <- unlist(list_results_all_chunks, recursive = FALSE)

    results_list[[i]] <- data.table::rbindlist(
      list_results_all_fragments,
      use.names = TRUE,
      fill = TRUE
    )
  }

  df_fragments_info_final <- data.table::rbindlist(
    Filter(Negate(is.null), results_list),
    use.names = TRUE, fill = TRUE
  )

  if (nrow(df_fragments_info_final) == 0) {
    stop("The final fRagmentomic dataframe is empty.")
  }

  if (verbose) {
    message(sprintf(
      "Processed %d fragment(s) from %d region(s).",
      nrow(df_fragments_info_final), nrow(df_regions)
    ))
  }

  if (!retain_fail_qc) {
    mask_fail_qc <- df_fragments_info_final[["Fragment_QC"]] != "OK"
    nfrag_fail_qc <- sum(mask_fail_qc)
    df_fragments_info_final <- df_fragments_info_final[!mask_fail_qc, , drop = FALSE]
    if (verbose) {
      message(sprintf(
        "Removed '%d' fragments that fail quality checks from output. Set `retain_fail_qc` to TRUE to retain them.",
        nfrag_fail_qc
      ))
    }
  }

  if (!is.null(output_path) && !is.na(output_path) && output_path != "") {
    if (file.exists(output_path)) {
      if (verbose) {
        message(sprintf("File '%s' already exists and will be overwritten.", output_path))
      }
    } else {
      output_parent <- dirname(output_path)
      if (!dir.exists(output_parent)) {
        if (verbose) {
          message(sprintf("Folder '%s' does not exist and will be created.", output_parent))
        }
        dir.create(output_parent, showWarnings = FALSE, recursive = TRUE)
      }
    }

    if (verbose) {
      message(sprintf("Writing results to: %s.", output_path))
    }
    write.table(df_fragments_info_final,
      file = output_path, sep = "\t", quote = FALSE,
      row.names = FALSE
    )

    return(NULL)
  }

  return(DataFrame(df_fragments_info_final))
}
