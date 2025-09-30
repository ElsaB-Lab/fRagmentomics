#' Plot Overall Nucleotide Frequency
#'
#' @description Creates a faceted bar plot to compare the overall proportion of each nucleotide across different groups. Includes
#' error bars (95% confidence intervals) and an optional global Chi-squared test for statistical significance.
#'
#' @param df_fragments The input dataframe containing fragment sequence data.
#' @param end_motif_5p Character string. Column name for 5' end sequences.
#' @param end_motif_3p Character string. Column name for 3' end sequences.
#' @param motif_type Character string. Which ends to analyze: 'Start', 'End', or 'Both'.
#' @param motif_size A single integer specifying the length of the motif to analyze.
#' @param col_z Character string. Column name for grouping. If NULL, no grouping is applied.
#' @param vals_z A character vector of group names from 'col_z' to include.
#'   If NULL, all unique groups in 'col_z' are used.
#' @param ... Additional arguments passed on to 'ggplot2::geom_bar()'.
#' @param colors_z A character vector of colors for the groups, or a single string
#'   naming an RColorBrewer palette (e.g., "Set2").
#' @param sample_id Sample identifier.
#' @param output_path Character vector for the plot output path.
#' @param ggsave_params A named list of arguments to be passed to 'ggplot2::ggsave()'. For example,
#'   'list(width = 8, height = 6, units = "in", dpi = 300, bg = "white")'. If not provided, sensible defaults will be used.
#'
#' @return A 'ggplot' object.
#'
#' @importFrom dplyr %>% filter select mutate all_of group_by summarise ungroup n bind_rows distinct
#' @importFrom stringr str_sub str_count
#' @importFrom tidyr pivot_longer
#' @importFrom stats prop.test chisq.test xtabs na.omit
#' @importFrom purrr map2 map_dbl
#' @importFrom scales percent_format pvalue
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' ## --- Create a dataset for demonstration ---
#' # Set a seed for reproducibility
#' set.seed(42)
#'
#' # Helper function to generate random DNA sequences with a bias
#' generate_biased_dna <- function(n_seq, len, prob) {
#'   bases <- c("A", "C", "G", "T")
#'   replicate(n_seq, paste(sample(bases, len, replace = TRUE, prob = prob), collapse = ""))
#' }
#'
#' # Create 50 "MUT" fragments with a high proportion of 'C'
#' df_mut <- data.frame(
#'   Fragment_Bases_5p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
#'   Fragment_Bases_3p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
#'   Fragment_Status_Simple = "MUT"
#' )
#'
#' # Create 50 "WT" fragments with a high proportion of 'G'
#' df_wt <- data.frame(
#'   Fragment_Bases_5p = generate_biased_dna(50, 10, prob = c(0.15, 0.15, 0.5, 0.2)),
#'   Fragment_Bases_3p = generate_biased_dna(50, 10, prob = c(0.15, 0.15, 0.5, 0.2)),
#'   Fragment_Status_Simple = "WT"
#' )
#'
#' # Combine into a single dataframe
#' example_df <- rbind(df_mut, df_wt)
#'
#' ## --- Function Calls ---
#'
#' # 1. Default plot: Compares MUT vs. WT groups for 3-mers
#' #    from both the 5' and 3' ends.
#' p1 <- plot_freq_barplot(example_df)
#' print(p1)
#'
#' # 2. Customized plot: Analyzes only the first nucleotide ('motif_size = 1')
#' #    of the 5' end ('motif_type = "Start"') using custom colors.
#' p2 <- plot_freq_barplot(
#'   df_fragments = example_df,
#'   motif_type = "Start",
#'   motif_size = 1,
#'   colors_z = c("MUT" = "#d95f02", "WT" = "#1b9e77")
#' )
#' print(p2)
#'
#' # 3. Ungrouped plot: Analyzes the overall nucleotide frequency
#' #    across all fragments combined.
#' p3 <- plot_freq_barplot(example_df, col_z = NULL)
#' print(p3)
#'
#' # 4. Plot with a subset of groups: If you had more than two groups
#' #    (e.g., "MUT", "WT", "AMB"), you could select specific ones to plot.
#' p4 <- plot_freq_barplot(
#'   df_fragments = example_df,
#'   vals_z = c("MUT", "WT")
#' )
#' print(p4)
#'
#' # 5. Save the default plot to a temporary folder.
#' # plot_freq_barplot(
#' #   df_fragments = example_df,
#' #   sample_id = "test01_freq",
#' #   output_path = "test01_freq_nucleotide_frequency.png"
#' # )
#'
#' # 6. Save a customized plot with specific dimensions.
#' # plot_freq_barplot(
#' #   df_fragments = example_df,
#' #   motif_type = "Start",
#' #   motif_size = 1,
#' #   sample_id = "test02_freq_custom",
#' #   output_path = "test02_freq_custom_nucleotide_frequency.png"
#' #   ggsave_params = list(width = 7, height = 5, units = "in")
#' # )
#'
plot_freq_barplot <- function(df_fragments,
                              end_motif_5p = "Fragment_Bases_5p",
                              end_motif_3p = "Fragment_Bases_3p",
                              motif_type = "Both",
                              motif_size = 3,
                              col_z = "Fragment_Status_Simple",
                              vals_z = NULL,
                              ...,
                              colors_z = "Set2",
                              sample_id = NA_character_,
                              output_path = NA_character_,
                              ggsave_params = list()) {
  # --- 1. Input Validation and Setup ---
  if (is.null(col_z) && !is.null(vals_z)) stop("If 'col_z' is NULL, 'vals_z' must also be NULL.")
  if (!is.null(col_z) && !col_z %in% names(df_fragments)) {
    stop(sprintf("Column '%s' not found in the dataframe.", col_z))
  }
  if (!motif_type %in% c("Start", "End", "Both")) stop("motif_type must be one of 'Start', 'End', or 'Both'.")
  is_grouped_analysis <- !is.null(col_z)
  if (!is_grouped_analysis) {
    df_fragments$placeholder_group <- "All Fragments"
    col_z <- "placeholder_group"
  }
  df_filtered <- df_fragments %>% filter(!is.na(.data[[col_z]]))
  if (is_grouped_analysis && !is.null(vals_z)) {
    df_filtered <- df_filtered %>% filter(.data[[col_z]] %in% vals_z)
  }
  if (nrow(df_filtered) == 0) stop("No data remains after filtering. Check 'col_z' and 'vals_z'.")
  if (is.null(vals_z)) {
    vals_z <- sort(unique(df_filtered[[col_z]]))
  }
  df_filtered[[col_z]] <- factor(df_filtered[[col_z]], levels = vals_z)
  if (is_grouped_analysis && length(vals_z) < 2) {
    stop("Grouped analysis requires at least two groups for comparison. Use col_z = NULL for a single group.")
  }

  # --- 2. Motif Size Adjustment ---
  max_len <- Inf
  if (motif_type %in% c("Start", "Both")) {
    max_len <- min(max_len, min(nchar(stats::na.omit(df_filtered[[end_motif_5p]]))))
  }
  if (motif_type %in% c("End", "Both")) {
    max_len <- min(max_len, min(nchar(stats::na.omit(df_filtered[[end_motif_3p]]))))
  }
  if (motif_size > max_len) {
    warning(sprintf(
      "Requested 'motif_size' (%d) is larger than the shortest available sequence (%d). Using maximum possible size: %d.",
      motif_size,
      max_len,
      max_len
    ))
    motif_size <- max_len
  }

  # --- 3. Data Preparation ---
  motifs_list <- list()
  if (motif_type %in% c("Start", "Both")) {
    motifs_list$start <- df_filtered %>%
      select(group = all_of(col_z), motif = all_of(end_motif_5p)) %>%
      mutate(motif = str_sub(motif, 1, motif_size))
  }
  if (motif_type %in% c("End", "Both")) {
    motifs_list$end <- df_filtered %>%
      select(group = all_of(col_z), motif = all_of(end_motif_3p)) %>%
      mutate(motif = str_sub(motif, -motif_size, -1))
  }
  nucleotide_cols <- c("A", "C", "G", "T")

  count_df <- bind_rows(motifs_list) %>%
    filter(!is.na(motif)) %>%
    mutate(
      A = str_count(motif, "A"), C = str_count(motif, "C"),
      G = str_count(motif, "G"), T = str_count(motif, "T")
    ) %>%
    select(-motif) %>%
    pivot_longer(
      cols = all_of(nucleotide_cols),
      names_to = "nucleotide",
      values_to = "count"
    ) %>%
    group_by(group, nucleotide) %>%
    summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop")

  # --- 4. Calculate Frequencies and Confidence Intervals ---
  plot_data <- count_df %>%
    group_by(group) %>%
    mutate(total_bases_in_group = sum(total_count)) %>%
    ungroup() %>%
    mutate(
      frequency = total_count / total_bases_in_group,
      ci = map2(total_count, total_bases_in_group, ~ stats::prop.test(.x, .y)$conf.int),
      ci_low = map_dbl(ci, ~ .x[1]),
      ci_high = map_dbl(ci, ~ .x[2])
    ) %>%
    mutate(
      nucleotide = factor(nucleotide, levels = c("A", "C", "G", "T")),
      group = factor(group, levels = vals_z)
    )

  # --- 5. Statistical Test & Plot Customization ---

  # Calculate total N for each GROUP and create new labels for the legend/axis
  group_totals <- plot_data %>%
    distinct(group, total_bases_in_group)

  new_group_labels <- paste0(group_totals$group, " (N=", group_totals$total_bases_in_group, ")")

  # Perform Chi-squared test
  plot_caption <- "Bars represent overall proportion with 95% confidence intervals."
  if (is_grouped_analysis) {
    contingency_table <- xtabs(total_count ~ nucleotide + group, data = plot_data)
    chi_test_result <- stats::chisq.test(contingency_table, simulate.p.value = TRUE)
    plot_caption <- paste0(plot_caption, "\nGlobal comparison (Chi-squared test), p-value = ", scales::pvalue(chi_test_result$p.value))
  }

  # Prepare color palette
  if (length(colors_z) == 1 && colors_z %in% row.names(RColorBrewer::brewer.pal.info)) {
    colors_z <- RColorBrewer::brewer.pal(n = length(vals_z), name = colors_z)
  }

  # --- 6. Generate Plot ---
  final_plot <- ggplot(plot_data, aes(x = group, y = frequency, fill = group)) +
    geom_bar(stat = "identity", width = 0.7, ...) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.25) +
    facet_wrap(~nucleotide, scales = "free_x", nrow = 1) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +

    # Use the new group labels in the scales for axis and legend
    scale_fill_manual(values = colors_z, name = "Group", labels = new_group_labels) +
    theme_bw(base_size = 14) +
    labs(
      title = paste0("Overall Nucleotide Frequency (", motif_size, "-mers)"),
      x = NULL,
      y = "Overall Frequency",
      caption = plot_caption
    ) +
    theme(
      plot.caption = element_text(hjust = 0, size = 10, face = "italic"),
      strip.text = element_text(face = "bold", size = 12),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_line(),
      legend.position = "right"
    )

  # --- 7. Save the plot to a file if an output_path is provided ---
  if (is.null(output_path) ||
      (is.character(output_path) && length(output_path) == 1 &&
      (is.na(output_path) || !nzchar(output_path)))) {
    return(final_plot)
  }

  if (!(is.character(output_path) && length(output_path) == 1)) {
    stop("'output_path' must be a single character string.")
  }

  save_plot_if_needed(final_plot, output_path, ggsave_params)
  return(invisible(NULL))
}
